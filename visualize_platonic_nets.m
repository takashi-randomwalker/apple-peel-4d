(* ================================================================
   visualize_platonic_nets.m
   Apple-Peel 展開図の可視化

   展開アルゴリズム:
   1. F1の法線が+z軸と平行になるよう全頂点を回転 (F1がxy平面に乗る)
   2. k=2..n: F(k-1)とFkの共有辺を軸にFkをxy平面に倒す
      倒す向きはF(k-1)の内側と反対側
   3. z座標を除いて2D座標を返す

   使い方:
     Get["/path/to/visualize_platonic_nets.m"]
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
If[scriptDir === "", scriptDir = NotebookDirectory[]];
Get[scriptDir <> "peeling3DLoxo.m"];
Get[scriptDir <> "dataPlatonic.mx"];

(* ----------------------------------------------------------------
   補助: 単位ベクトル from を to に重ねる回転行列 (Rodrigues の公式)
   ---------------------------------------------------------------- *)
rotationMatrix3D[from_, to_] :=
  Module[{v, c, s, vx, perp},
    v = Cross[from, to];
    c = from . to;
    s = Norm[v];
    Which[
      s < 10^-10 && c > 0,   (* 同方向 *)
        IdentityMatrix[3],
      s < 10^-10 && c <= 0,  (* 反対方向: from に垂直な軸で 180 度回転 *)
        perp = If[Abs[from[[1]]] < 0.9,
                  Normalize[Cross[from, {1., 0., 0.}]],
                  Normalize[Cross[from, {0., 1., 0.}]]];
        2 Outer[Times, perp, perp] - IdentityMatrix[3],
      True,
        vx = {{0, -v[[3]], v[[2]]},
              {v[[3]], 0, -v[[1]]},
              {-v[[2]], v[[1]], 0}};
        IdentityMatrix[3] + vx + vx . vx (1 - c) / s^2
    ]
  ];

(* ----------------------------------------------------------------
   2D展開図生成: p3lUnfoldNet

   各面に独立した2D座標辞書を持つ．同じ頂点でも面ごとに異なる
   2D位置を持てるため，展開図が正しく生成される．

   pv0 = F1 を xy 平面に乗せた後の3D座標（固定参照）
   facePos[[k]] = face order[[k]] の 頂点 -> 2D座標 の辞書
   ---------------------------------------------------------------- *)
p3lUnfoldNet[vers_, faces_, order_] :=
  Module[{nf, pv0, verForFaces, normals, changed,
          baseFaceNum, prevFaceNum, n1, n2, axis, angle,
          edgeVerts, pivotPt, rotFn, lastNormal, R,
          polys2D, dir, a},
    nf          = Length[order];
    pv0         = N[vers];
    verForFaces = Map[pv0[[#]] &, faces, {2}];

    (* 各面の外向き単位法線 *)
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

      (* F_{k-1} を元の 3D 位置のまま追加 *)
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


(* ----------------------------------------------------------------
   グラデーション色: 青(先頭) → 緑(中間) → 赤(末尾)
   ---------------------------------------------------------------- *)
peelingColor[k_, nfTotal_] :=
  Blend[{RGBColor[0.20, 0.45, 0.80],
         RGBColor[0.15, 0.70, 0.40],
         RGBColor[0.90, 0.25, 0.15]},
        (k - 1) / Max[nfTotal - 1, 1]];


(* ----------------------------------------------------------------
   成功展開図の描画
   ---------------------------------------------------------------- *)
netGraphicsOK[polys2D_, nfTotal_, sz_: 140] :=
  Graphics[
    Table[{FaceForm[peelingColor[k, nfTotal]],
           EdgeForm[Directive[GrayLevel[0.3], Thickness[0.005]]],
           Polygon[polys2D[[k]]]},
          {k, Length[polys2D]}],
    ImageSize -> sz,
    PlotRangePadding -> Scaled[0.05],
    Background -> GrayLevel[0.97]];


(* ----------------------------------------------------------------
   失敗展開図の描画 (ピンク背景・赤枠・× マーク)
   ---------------------------------------------------------------- *)
netGraphicsFail[polys2D_, nPlaced_, nfTotal_, sz_: 140] :=
  Module[{allPts, xMin, xMax, yMin, yMax, xMid, yMid, r, pad},
    allPts = Flatten[polys2D, 1];
    xMin = Min[allPts[[All, 1]]]; xMax = Max[allPts[[All, 1]]];
    yMin = Min[allPts[[All, 2]]]; yMax = Max[allPts[[All, 2]]];
    xMid = (xMin + xMax) / 2;
    yMid = (yMin + yMax) / 2;
    r    = Max[xMax - xMin, yMax - yMin] / 2;
    pad  = 0.30 r;  (* テキスト・×マーク用の余白 *)
    Graphics[
      {Table[{FaceForm[peelingColor[k, nfTotal]],
              EdgeForm[Directive[GrayLevel[0.4], Thickness[0.005]]],
              Polygon[polys2D[[k]]]},
             {k, Length[polys2D]}],
       {Thickness[0.025], RGBColor[0.85, 0.10, 0.10],
        Line[{{xMid + 0.55*r, yMid + 0.55*r}, {xMid + 0.85*r, yMid + 0.85*r}}],
        Line[{{xMid + 0.85*r, yMid + 0.55*r}, {xMid + 0.55*r, yMid + 0.85*r}}]},
       Text[Style[ToString[nPlaced] <> "/" <> ToString[nfTotal],
                  FontSize -> 10, FontColor -> RGBColor[0.7, 0., 0.], FontWeight -> Bold],
            {xMid, yMin - pad}]},
      ImageSize -> sz,
      PlotRange -> {{xMin - pad, xMax + pad}, {yMin - 2.5*pad, yMax + pad}},
      Background -> GrayLevel[0.97],
      PlotRangePadding -> None]
  ];


(* ----------------------------------------------------------------
   メインループ: 全出力を Column にまとめて返す
   ---------------------------------------------------------------- *)
solidNames   = {"Tetrahedron", "Cube", "Octahedron",
                "Dodecahedron", "Icosahedron"};
ruleKeys     = {"R1", "R2", "R3"};
ruleLabels   = {"R1 (max \[Phi])", "R2 (loxodrome)", "R3 (max z)"};
fbKeys       = {"withFallback", "noFallback"};
fbLabels     = {"with fallback", "no fallback"};
netSize      = 120;
maxShowOK    = Infinity;
maxShowFail  = Infinity;
(* True: 同じペーリング順序 (order) を持つ結果を1件に絞る
   False: 全ペア (top, f2) の結果をすべて表示 *)
deduplicate  = True;

(* 1ペア分の成功/失敗行を output に追記するヘルパー *)
appendNetRows[vers_, faces_, nf_, nTotal_, allRes_, label_] :=
  Module[{complete, failed, uniqOK, uniqFail, sampleOK, sampleFail, nComp},
    complete = Select[allRes,  #[[4]] &];
    failed   = Select[allRes, !#[[4]] &];
    nComp    = Length[complete];

    (* deduplicate = True のとき: 同じ order を持つ結果を1件に絞る *)
    uniqOK   = If[deduplicate, DeleteDuplicatesBy[complete, #[[3]] &], complete];
    uniqFail = If[deduplicate, DeleteDuplicatesBy[failed,   #[[3]] &], failed];

    sampleOK   = Take[uniqOK,   UpTo[maxShowOK]];
    sampleFail = Take[uniqFail, UpTo[maxShowFail]];

    AppendTo[output,
      Style["    [" <> label <> "]  \:6210\:529f: " <> ToString[nComp] <>
            "/" <> ToString[nTotal] <>
            "   \:5931\:6557: " <> ToString[nTotal - nComp] <> "/" <> ToString[nTotal],
            12, RGBColor[0.35, 0.35, 0.55]]];

    If[Length[sampleOK] > 0,
      AppendTo[output,
        Row[(netGraphicsOK[p3lUnfoldNet[vers, faces, #[[3]]], nf, netSize]) & /@
            sampleOK, Spacer[3]]]
    ];

    If[Length[sampleFail] > 0,
      AppendTo[output,
        Row[(netGraphicsFail[p3lUnfoldNet[vers, faces, #[[3]]],
                             Length[#[[3]]], nf, netSize]) & /@
            sampleFail, Spacer[3]]]
    ]
  ];

output = {};

Do[
  {vers, faces} = p3lExtractPolyhedron[name];
  nf     = platonicData[name]["nFaces"];
  nTotal = platonicData[name]["totalPairs"];

  AppendTo[output, Style[StringRepeat["\[LongDash]", 50], Gray]];
  AppendTo[output,
    Style[name <> "  (faces: " <> ToString[nf] <> ")", Bold, 16]];

  Do[
    AppendTo[output,
      Style["  " <> ruleLabels[[ri]], 14, RGBColor[0.2, 0.2, 0.6]]];
    Do[
      appendNetRows[vers, faces, nf, nTotal,
        platonicData[name][ruleKeys[[ri]]][fbKeys[[fi]]],
        fbLabels[[fi]]],
      {fi, 2}
    ],
    {ri, 3}
  ],
  {name, solidNames}
];

(* カラーバー *)
AppendTo[output, Style[StringRepeat["\[LongDash]", 50], Gray]];
AppendTo[output,
  Style["\:30ab\:30e9\:30fc\:30d0\:30fc: \:9752=\:5148\:982d\:306e\:9762  \:7dd1=\:4e2d\:9593  \:8d64=\:672b\:5c3e\:306e\:9762" <>
        "  (\:5931\:6557\:306f\:30d4\:30f3\:30af\:80cc\:666f\:30fb\:8d64\:679d\[Times])",
        11, Gray]];
AppendTo[output,
  Row[Table[
    Graphics[{FaceForm[peelingColor[k, 20]], EdgeForm[None],
              Rectangle[{k - 1, 0}, {k, 1}]},
             ImageSize -> {300, 16},
             PlotRange -> {{0, 20}, {0, 1}},
             AspectRatio -> Full],
    {k, 1, 20}],
  Spacer[0]]];

Column[output, Spacings -> 1]
