(* ================================================================
   gen_exampleOfPeeling_v2.m
   切頂二十面体 (TruncatedIcosahedron) の Apple-Peel スナップショット図
   Method B: 蓄積面を F_k の平面に倒す (正確な二面角使用)
   19面目スナップショット，45度傾け，赤/緑/黄/灰色
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
If[scriptDir === "", scriptDir = NotebookDirectory[]];
Get[scriptDir <> "peeling3DLoxo.m"];

(* ----------------------------------------------------------------
   補助: 回転行列 (Rodrigues)  ※ visualize_*.m から移植
   ---------------------------------------------------------------- *)
rotationMatrix3D[from_, to_] :=
  Module[{v, c, s, vx, perp},
    v = Cross[from, to]; c = from . to; s = Norm[v];
    Which[
      s < 10^-10 && c > 0, IdentityMatrix[3],
      s < 10^-10 && c <= 0,
        perp = If[Abs[from[[1]]] < 0.9,
                  Normalize[Cross[from, {1., 0., 0.}]],
                  Normalize[Cross[from, {0., 1., 0.}]]];
        2 Outer[Times, perp, perp] - IdentityMatrix[3],
      True,
        vx = {{0, -v[[3]], v[[2]]}, {v[[3]], 0, -v[[1]]}, {-v[[2]], v[[1]], 0}};
        IdentityMatrix[3] + vx + vx . vx (1 - c) / s^2
    ]];

(* ----------------------------------------------------------------
   2D net 生成: Method B (蓄積面を F_k の平面に倒す，正確な二面角使用)
   ---------------------------------------------------------------- *)
p3lUnfoldNet[vers_, faces_, order_] :=
  Module[{nf, pv0, verForFaces, normals, changed,
          baseFaceNum, prevFaceNum, n1, n2, axis, angle,
          edgeVerts, pivotPt, rotFn, lastNormal, R,
          polys2D, dir, a},
    nf          = Length[order];
    pv0         = N[vers];
    verForFaces = Map[pv0[[#]] &, faces, {2}];
    normals = Table[
      With[{verts = verForFaces[[i]]},
        With[{nn = Normalize[Cross[verts[[2]] - verts[[1]],
                                   verts[[3]] - verts[[1]]]]},
          If[nn . Mean[verts] < 0, -nn, nn]]],
      {i, Length[faces]}];
    changed = {};
    Do[
      baseFaceNum = order[[k]];
      prevFaceNum = order[[k - 1]];
      AppendTo[changed, verForFaces[[prevFaceNum]]];
      n1    = normals[[prevFaceNum]];
      n2    = normals[[baseFaceNum]];
      axis  = Cross[n1, n2];
      If[Norm[axis] < 10^-10, Continue[]];
      axis  = Normalize[axis];
      angle = ArcCos[Clip[n1 . n2, {-1, 1}]];
      edgeVerts = Intersection[faces[[prevFaceNum]], faces[[baseFaceNum]]];
      pivotPt   = pv0[[ edgeVerts[[1]] ]];
      rotFn   = RotationTransform[angle, axis, pivotPt];
      changed = Map[rotFn, changed, {2}],
      {k, 2, nf}
    ];
    AppendTo[changed, verForFaces[[ order[[-1]] ]]];
    lastNormal = normals[[ order[[-1]] ]];
    R          = rotationMatrix3D[lastNormal, {0., 0., 1.}];
    polys2D    = Map[((# . Transpose[R])[[{1, 2}]] &), changed, {2}];
    dir = Mean[polys2D[[2]]] - Mean[polys2D[[1]]];
    a   = ArcTan[dir[[1]], dir[[2]]];
    Map[Function[pt, {pt[[1]] Cos[a] + pt[[2]] Sin[a],
                     -pt[[1]] Sin[a] + pt[[2]] Cos[a]}],
        polys2D, {2}]
  ];

(* ================================================================
   データ読み込みとペーリング
   ================================================================ *)
{vers, faces} = p3lExtractPolyhedron["TruncatedIcosahedron"];
nf  = Length[faces];
top = 1;

(* Face-up 座標 (3D 可視化用) *)
lv = p3lAlignTopToZ[vers, faces, top];
cc = Table[p3lFaceCentroid[i, lv, faces], {i, nf}];

(* RZ ペーリング: 全ペア計算 *)
Print["Computing RZ peeling for TruncatedIcosahedron..."];
results = p3lRunAll[vers, faces, "maxz"];

(* top=1 の最初の成功ペアの order を取得 *)
r1pairs = Select[results, #[[1]] == top && #[[4]] == True &];
If[Length[r1pairs] == 0, Print["ERROR: no successful pair for top=1"]; Abort[]];
order = r1pairs[[1, 3]];
Print["Full order length: ", Length[order], "  (expected 32)"];

(* 19面目スナップショット *)
snapStep = 19;
orderSnap = order[[1 ;; snapStep]];

(* ================================================================
   2D Net (灰色)
   ================================================================ *)
polys2D = p3lUnfoldNet[vers, faces, orderSnap];

netGraphic = Graphics[
  {GrayLevel[0.55],
   EdgeForm[{GrayLevel[0.2], Thickness[0.003]}],
   Polygon /@ polys2D},
  Background -> White,
  PlotRangePadding -> Scaled[0.08]
];

(* ================================================================
   3D 多面体
   ================================================================ *)

(* 赤面: ペーリング順の最終面 (RZ では対蹠面と一致) *)
lastFace  = order[[-1]];
peeledSet = AssociationThread[orderSnap -> Range[snapStep]];

verForFaces = Map[lv[[#]] &, faces, {2}];

facePolygons = Table[
  col = Which[
    i == lastFace,            RGBColor[0.90, 0.15, 0.15],   (* 赤: 最終面 *)
    KeyExistsQ[peeledSet, i], RGBColor[0.88, 0.76, 0.10],   (* 黄: 剥き済み *)
    True,                     RGBColor[0.20, 0.72, 0.25]    (* 緑: 未剥き *)
  ];
  {col, EdgeForm[{GrayLevel[0.25], Thickness[0.001]}],
   Polygon[verForFaces[[i]]]},
  {i, nf}];

(* c_1 方向の矢印: 重心から c_1 方向（+z）の外側へ向かう
   旧コード: Arrow[{Mean[vers], 1.7 * fCenters[[order[[1]]]]}] に対応 *)
c1Arrow = {Arrowheads[0.05], AbsoluteThickness[2], Black,
           Arrow[{Mean[lv], 1.7 cc[[top]]}]};

(* 45度傾いた視点: 仰角45度，やや右前方から *)
(* z成分とxy成分を等しく取ると仰角arctan(z/sqrt(x^2+y^2))=45度 *)
poly3D = Graphics3D[
  {facePolygons, c1Arrow},
  ViewPoint  -> {2.0, -2.0, 2.0},
  ViewVertical -> {0, 0, 1},
  Boxed      -> False,
  Lighting   -> "Neutral",
  SphericalRegion -> True
];

(* ================================================================
   組み合わせて PDF 出力
   ================================================================ *)
combined = Row[{netGraphic, Spacer[20], poly3D},
               BaselinePosition -> Center];

outFile = scriptDir <> "exampleOfPeeling_v2.pdf";
Export[outFile, combined, ImageResolution -> 150];
Print["Exported: ", outFile];
