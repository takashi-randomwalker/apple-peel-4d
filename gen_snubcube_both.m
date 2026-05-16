(* ================================================================
   gen_snubcube_both.m
   Snub Cube オリジナルと鏡像の RZ 展開図比較図を生成
   - オリジナル: {snub,snub} + {gyrate,snub} 軌道が成功
   - 鏡像:       {snub,snub} + {sq,snub} 軌道が成功（正方形スタート）
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
If[scriptDir === "", scriptDir = NotebookDirectory[]];
Get[scriptDir <> "peeling3DLoxo.m"];

(* ----------------------------------------------------------------
   回転行列 (Rodrigues)
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
   2D 展開図生成 (Method B)
   ---------------------------------------------------------------- *)
p3lUnfoldNet[vers_, faces_, order_] :=
  Module[{nf, pv0, verForFaces, normals, changed,
          baseFaceNum, prevFaceNum, n1, n2, axis, angle,
          edgeVerts, pivotPt, rotFn, lastNormal, R,
          polys2D, dir, a},
    nf = Length[order]; pv0 = N[vers];
    verForFaces = Map[pv0[[#]] &, faces, {2}];
    normals = Table[
      With[{verts = verForFaces[[i]]},
        With[{nn = Normalize[Cross[verts[[2]] - verts[[1]],
                                   verts[[3]] - verts[[1]]]]},
          If[nn . Mean[verts] < 0, -nn, nn]]],
      {i, Length[faces]}];
    changed = {};
    Do[
      baseFaceNum = order[[k]]; prevFaceNum = order[[k - 1]];
      AppendTo[changed, verForFaces[[prevFaceNum]]];
      n1 = normals[[prevFaceNum]]; n2 = normals[[baseFaceNum]];
      axis = Cross[n1, n2];
      If[Norm[axis] < 10^-10, Continue[]];
      axis = Normalize[axis];
      angle = ArcCos[Clip[n1 . n2, {-1, 1}]];
      edgeVerts = Intersection[faces[[prevFaceNum]], faces[[baseFaceNum]]];
      pivotPt = pv0[[edgeVerts[[1]]]];
      rotFn = RotationTransform[angle, axis, pivotPt];
      changed = Map[rotFn, changed, {2}],
      {k, 2, nf}
    ];
    AppendTo[changed, verForFaces[[order[[-1]]]]];
    lastNormal = normals[[order[[-1]]]];
    R = rotationMatrix3D[lastNormal, {0., 0., 1.}];
    polys2D = Map[((# . Transpose[R])[[{1, 2}]] &), changed, {2}];
    dir = Mean[polys2D[[2]]] - Mean[polys2D[[1]]];
    a = ArcTan[dir[[1]], dir[[2]]];
    Map[Function[pt, {pt[[1]] Cos[a] + pt[[2]] Sin[a],
                     -pt[[1]] Sin[a] + pt[[2]] Cos[a]}],
        polys2D, {2}]
  ];

(* ----------------------------------------------------------------
   描画ヘルパー
   ---------------------------------------------------------------- *)
peelingColor[k_, n_] := GrayLevel[0.2 + 0.5 (k - 1) / Max[n - 1, 1]];

makeNetGraphic[vers_, faces_, order_, sz_: 150] :=
  Module[{polys2D = p3lUnfoldNet[vers, faces, order], nTotal = Length[order]},
    Graphics[
      Table[{FaceForm[peelingColor[k, nTotal]],
             EdgeForm[Directive[Black, Thickness[0.004]]],
             Polygon[polys2D[[k]]]},
            {k, nTotal}],
      ImageSize -> {sz, sz},
      PlotRangePadding -> Scaled[0.08],
      Background -> White]
  ];

(* ----------------------------------------------------------------
   Snub Cube データと面タイプ分類
   ---------------------------------------------------------------- *)
{vers, faces} = p3lExtractPolyhedron["SnubCube"];
nf       = Length[faces];
adj      = p3lBuildAdj[faces];
faceType = Length /@ faces;

squareFaces = Select[Range[nf], faceType[[#]] == 4 &];
triFaces    = Select[Range[nf], faceType[[#]] == 3 &];

(* snub 三角形: 正方形と1辺共有，gyrate 三角形: 正方形と非隣接 *)
triAdjSqCount = Association @ Table[
  f -> Length[Intersection[adj[[f]], squareFaces]], {f, triFaces}];
snubTri   = Select[triFaces, triAdjSqCount[#] == 1 &];
gyrateTri = Select[triFaces, triAdjSqCount[#] == 0 &];

getSubtype[f_] := Which[
  faceType[[f]] == 4,   "sq",
  MemberQ[snubTri, f],  "snub",
  True,                 "gyrate"
];

Print["Face types — square: ", Length[squareFaces],
      "  snub-tri: ", Length[snubTri],
      "  gyrate-tri: ", Length[gyrateTri]];

(* ----------------------------------------------------------------
   RZ 計算: オリジナルと鏡像
   ---------------------------------------------------------------- *)
versMirror = vers * ConstantArray[{-1., 1., 1.}, Length[vers]];

Print["Computing RZ for original Snub Cube..."];
r3orig = p3lRunAll[vers,       faces, "maxz"];
Print["Computing RZ for mirror Snub Cube..."];
r3mir  = p3lRunAll[versMirror, faces, "maxz"];

succOrig = Select[r3orig, #[[4]] &];
succMir  = Select[r3mir,  #[[4]] &];

Print["Original success: ", Length[succOrig], "/", Length[r3orig]];
Print["Mirror  success: ", Length[succMir],  "/", Length[r3mir]];

(* サブ軌道ごとの集計 *)
tally[results_] := Tally[{getSubtype[#[[3,1]]], getSubtype[#[[3,2]]]} & /@
                          Select[results, #[[4]] &]];
Print["Original sub-orbits: ", tally[r3orig]];
Print["Mirror   sub-orbits: ", tally[r3mir]];

(* ----------------------------------------------------------------
   各サブ軌道の代表を選択（最初の成功例1件ずつ）
   ---------------------------------------------------------------- *)
pickRep[results_, f1type_, f2type_] :=
  Module[{filtered = Select[results,
            TrueQ[#[[4]]] &&
            getSubtype[#[[3, 1]]] === f1type &&
            getSubtype[#[[3, 2]]] === f2type &]},
    If[Length[filtered] == 0, None, filtered[[1]]]
  ];

repOrigSS = pickRep[r3orig, "snub",   "snub"];   (* {snub,snub}   *)
repOrigGS = pickRep[r3orig, "gyrate", "snub"];   (* {gyrate,snub} *)
repMirrSS = pickRep[r3mir,  "snub",   "snub"];   (* {snub,snub}   *)
repMirrSQ = pickRep[r3mir,  "sq",     "snub"];   (* {sq,snub}     *)

Print["Reps — orig snub-snub: ", repOrigSS =!= None,
      "  orig gyrate-snub: ", repOrigGS =!= None,
      "  mirr snub-snub: ", repMirrSS =!= None,
      "  mirr sq-snub: ", repMirrSQ =!= None];

(* ----------------------------------------------------------------
   図の組み立て
   ---------------------------------------------------------------- *)
blankCell = Graphics[{}, ImageSize -> {150, 180}, Background -> White];

netOrSS = If[ListQ[repOrigSS] && Length[repOrigSS] >= 3,
  makeNetGraphic[vers,       faces, repOrigSS[[3]], 150], blankCell];
netOrGS = If[ListQ[repOrigGS] && Length[repOrigGS] >= 3,
  makeNetGraphic[vers,       faces, repOrigGS[[3]], 150], blankCell];
netMiSS = If[ListQ[repMirrSS] && Length[repMirrSS] >= 3,
  makeNetGraphic[versMirror, faces, repMirrSS[[3]], 150], blankCell];
netMiSQ = If[ListQ[repMirrSQ] && Length[repMirrSQ] >= 3,
  makeNetGraphic[versMirror, faces, repMirrSQ[[3]], 150], blankCell];

fig = Grid[{
  {Style["", 1],
   Style["Original", Bold, 12, RGBColor[0.15, 0.25, 0.65]],
   Style["Original", Bold, 12, RGBColor[0.15, 0.25, 0.65]],
   Style["Mirror",   Bold, 12, RGBColor[0.15, 0.55, 0.25]],
   Style["Mirror",   Bold, 12, RGBColor[0.15, 0.55, 0.25]]},
  {Style["orbit", 9, Gray],
   Style[Row[{Subscript["F",1], " = snub \[CapitalDelta]"}],   9, Gray],
   Style[Row[{Subscript["F",1], " = gyrate \[CapitalDelta]"}], 9, Gray],
   Style[Row[{Subscript["F",1], " = snub \[CapitalDelta]"}],   9, Gray],
   Style[Row[{Subscript["F",1], " = \[FilledSquare]"}],         9, Gray]},
  {Style["RZ net", 9, Gray], netOrSS, netOrGS, netMiSS, netMiSQ}
  },
  Spacings -> {1, 0.5},
  Alignment -> Center,
  Background -> White
];

outFile = scriptDir <> "snubcube_both.pdf";
Export[outFile, Rasterize[fig, ImageResolution -> 250, Background -> White]];
Print["Exported: ", outFile];
