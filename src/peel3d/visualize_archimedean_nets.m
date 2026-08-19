(* ================================================================
   visualize_archimedean_nets.m
   Apple-Peel 展開図の可視化 — Archimedean Solids 版

   データソース:
     archimedean_faceup_results.mx
       allResults[name]["details"]["R1"|"R2"|"R3"]
         = {{top, f2, order, success}, ...}   (with fallback のみ)
     allResults[name]["nFaces"]      ... 面数
     allResults[name]["totalPairs"]  ... ペア総数

   表示規則:
     RS (Spiral rule, max phi, key "R1")
     RZ (Zonal rule,  max z,   key "R3")
     R2 は showR2 = True のときのみ表示

   使い方 (Notebook から):
     Get["/path/to/visualize_archimedean_nets.m"]
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
If[scriptDir === "", scriptDir = NotebookDirectory[]];
Get[scriptDir <> "peeling3DLoxo.m"];
Get[scriptDir <> "archimedean_faceup_results.mx"];
If[TrueQ[showNoFallback],
  Get[scriptDir <> "archimedean_nofallback_results.mx"]];

(* ----------------------------------------------------------------
   設定パラメータ（Get 前に設定済みの場合は上書きしない）
   ---------------------------------------------------------------- *)
If[!ValueQ[netSize],        netSize        = 100];
If[!ValueQ[maxShowOK],      maxShowOK      = Infinity];
If[!ValueQ[maxShowFail],    maxShowFail    = Infinity];
If[!ValueQ[deduplicate],    deduplicate    = True];
If[!ValueQ[deduplicateGeo], deduplicateGeo = False];
If[!ValueQ[showR2],          showR2          = False];
If[!ValueQ[showNoFallback],  showNoFallback  = False];

(* ----------------------------------------------------------------
   幾何学的重複除去キー
   各面の頂点数の列（ペーリング順）をキーとする．
   2D 座標計算が不要で高速．
   同じ多面体内の等変な軌道は同一キーになるため，
   面タイプが異なる軌道（三角形スタート vs 正方形スタート）を正しく区別する．
   ---------------------------------------------------------------- *)
netGeoKey[vers_, faces_, order_] :=
  Length[faces[[#]]] & /@ order;

(* ----------------------------------------------------------------
   補助: 単位ベクトル from を to に重ねる回転行列 (Rodrigues の公式)
   ---------------------------------------------------------------- *)
rotationMatrix3D[from_, to_] :=
  Module[{v, c, s, vx, perp},
    v = Cross[from, to];
    c = from . to;
    s = Norm[v];
    Which[
      s < 10^-10 && c > 0,
        IdentityMatrix[3],
      s < 10^-10 && c <= 0,
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
   2D 展開図生成: p3lUnfoldNet  (Method B — 厳密二面角)
   各ステップで F_1..F_{k-1} を F_k の平面に折り込み，
   最終的に F_n の平面に全面が収まった状態を 2D 射影する．
   切断不要: 既存の net を共有辺周りに転がす連続操作．
   回転角 = 外面二面角 = ArcCos[n_{k-1} · n_k]（厳密）
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
      axis  = Cross[n1, n2];          (* 共有辺方向に平行 *)
      If[Norm[axis] < 10^-10, Continue[]];
      axis  = Normalize[axis];
      angle = ArcCos[Clip[n1 . n2, {-1, 1}]];  (* 外面二面角（厳密）*)

      (* ピボット: 共有辺の頂点（元の 3D 位置） *)
      edgeVerts = Intersection[faces[[prevFaceNum]], faces[[baseFaceNum]]];
      pivotPt   = pv0[[ edgeVerts[[1]] ]];

      rotFn   = RotationTransform[angle, axis, pivotPt];
      changed = Map[rotFn, changed, {2}],
      {k, 2, nf}
    ];

    (* F_n を元の 3D 位置で追加（出力平面を定義） *)
    AppendTo[changed, verForFaces[[ order[[-1]] ]]];

    (* F_n の外向き法線を +z に揃えて 2D 射影 *)
    lastNormal = normals[[ order[[-1]] ]];
    R          = rotationMatrix3D[lastNormal, {0., 0., 1.}];
    polys2D    = Map[((# . Transpose[R])[[{1, 2}]] &), changed, {2}];

    (* F_1→F_2 方向を +x に揃える *)
    dir = Mean[polys2D[[2]]] - Mean[polys2D[[1]]];
    a   = ArcTan[dir[[1]], dir[[2]]];
    Map[Function[pt, {pt[[1]] Cos[a] + pt[[2]] Sin[a],
                     -pt[[1]] Sin[a] + pt[[2]] Cos[a]}],
        polys2D, {2}]
  ];

(* ----------------------------------------------------------------
   グラデーション色: 薄灰(先頭, 0.7) -> 濃灰(末尾, 0.2)
   ---------------------------------------------------------------- *)
peelingColor[k_, nfTotal_] :=
  GrayLevel[0.5 + 0.4 * (k - 1) / Max[nfTotal - 1, 1]];

(* ----------------------------------------------------------------
   成功展開図の描画
   ---------------------------------------------------------------- *)
netGraphicsOK[polys2D_, nfTotal_, sz_: 120] :=
  Graphics[
    Table[{FaceForm[peelingColor[k, nfTotal]],
           EdgeForm[Directive[Black, Thickness[0.005]]],
           Polygon[polys2D[[k]]]},
          {k, Length[polys2D]}],
    ImageSize -> sz,
    PlotRangePadding -> Scaled[0.05],
    Background -> GrayLevel[0.97]];

(* ----------------------------------------------------------------
   失敗展開図の描画 (赤 x マーク・配置数/全数ラベル)
   ---------------------------------------------------------------- *)
netGraphicsFail[polys2D_, nPlaced_, nfTotal_, sz_: 120] :=
  Module[{allPts, xMin, xMax, yMin, yMax, xMid, yMid, r, pad},
    allPts = Flatten[polys2D, 1];
    xMin = Min[allPts[[All, 1]]]; xMax = Max[allPts[[All, 1]]];
    yMin = Min[allPts[[All, 2]]]; yMax = Max[allPts[[All, 2]]];
    xMid = (xMin + xMax) / 2;
    yMid = (yMin + yMax) / 2;
    r    = Max[xMax - xMin, yMax - yMin] / 2;
    pad  = 0.30 r;
    Graphics[
      {Table[{FaceForm[peelingColor[k, nfTotal]],
              EdgeForm[Directive[Black, Thickness[0.005]]],
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
   1変種 (rule) の展開図行を output に追記するヘルパー
   allRes: {{top, f2, order, success}, ...}
   ---------------------------------------------------------------- *)
appendNetRowsArch[vers_, faces_, nf_, nTotal_, allRes_, label_] :=
  Module[{complete, failed, uniqOK, uniqFail, sampleOK, sampleFail, nComp},
    complete = Select[allRes,  #[[4]] &];
    failed   = Select[allRes, !#[[4]] &];
    nComp    = Length[complete];

    uniqOK   = If[deduplicate, DeleteDuplicatesBy[complete, #[[3]] &], complete];
    uniqFail = If[deduplicate, DeleteDuplicatesBy[failed,   #[[3]] &], failed];

    If[deduplicateGeo,
      uniqOK   = DeleteDuplicatesBy[uniqOK,   netGeoKey[vers, faces, #[[3]]] &];
      uniqFail = DeleteDuplicatesBy[uniqFail, netGeoKey[vers, faces, #[[3]]] &]
    ];

    sampleOK   = Take[uniqOK,   UpTo[maxShowOK]];
    sampleFail = Take[uniqFail, UpTo[maxShowFail]];

    geoNote = If[deduplicateGeo,
      "  (geo-uniq: " <> ToString[Length[uniqOK]] <> ")", ""];
    AppendTo[output,
      Style["    [" <> label <> "]  \:6210\:529f: " <> ToString[nComp] <>
            "/" <> ToString[nTotal] <>
            "   \:5931\:6557: " <> ToString[nTotal - nComp] <> "/" <> ToString[nTotal] <>
            geoNote,
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

(* ================================================================
   表示する規則のリスト
   ================================================================ *)
ruleKeyList       = {"R1", "R3"};
ruleShortNameList = {"RS", "RZ"};
ruleLabelList     = {"RS (Spiral rule, max \[Phi])", "RZ (Zonal rule, max z)"};
ruleColorList     = {RGBColor[0.2, 0.2, 0.6], RGBColor[0.1, 0.5, 0.1]};

If[showR2,
  ruleKeyList       = Append[ruleKeyList,       "R2"];
  ruleShortNameList = Append[ruleShortNameList,  "R2"];
  ruleLabelList     = Append[ruleLabelList,      "R2 (loxodrome, min \[Phi])"];
  ruleColorList     = Append[ruleColorList,      RGBColor[0.5, 0.2, 0.5]]
];

(* ================================================================
   アルキメデス13種のリスト
   ================================================================ *)
solidNames = {
  "TruncatedTetrahedron",
  "Cuboctahedron",
  "TruncatedCube",
  "TruncatedOctahedron",
  "Rhombicuboctahedron",
  "TruncatedCuboctahedron",
  "SnubCube",
  "Icosidodecahedron",
  "TruncatedDodecahedron",
  "TruncatedIcosahedron",
  "Rhombicosidodecahedron",
  "TruncatedIcosidodecahedron",
  "SnubDodecahedron"
};

(* ================================================================
   メインループ
   ================================================================ *)
output = {};

AppendTo[output,
  Style["Archimedean Solids — Apple-Peel Unfolding  (Face-up, with fallback)",
        Bold, 15, RGBColor[0.1, 0.3, 0.6]]];
AppendTo[output,
  Style["\:6a19\:6e96\:6761\:4ef6 (eps > 0, thresh = \[Minus]10^\[Minus]10)  \
\[Rule] det(c\:2081, c\:2096, c\:2097) \[GreaterEqual] \[Minus]\[Epsilon]",
        12, GrayLevel[0.4]]];

Do[
  {vers, faces} = p3lExtractPolyhedron[name];
  nf     = allResults[name]["nFaces"];
  nTotal = allResults[name]["totalPairs"];

  (* 成功率サマリ行を計算 *)
  summaryParts = Table[
    nSucc = allResults[name][ruleKeyList[[ri]]];
    pct   = NumberForm[N[100 nSucc / nTotal], {4, 1}];
    ruleShortNameList[[ri]] <> ": " <> ToString[nSucc] <> "/" <> ToString[nTotal] <>
    " (" <> ToString[pct] <> "%)",
    {ri, Length[ruleKeyList]}
  ];

  AppendTo[output, Style[StringRepeat["\[LongDash]", 55], Gray]];
  AppendTo[output,
    Style[name <> "  (faces: " <> ToString[nf] <> ",  pairs: " <> ToString[nTotal] <> ")",
          Bold, 16]];
  AppendTo[output,
    Style["    " <> StringRiffle[summaryParts, "   "], 12, GrayLevel[0.4]]];

  Do[
    AppendTo[output,
      Style["  " <> ruleLabelList[[ri]], 13, ruleColorList[[ri]]]];
    appendNetRowsArch[vers, faces, nf, nTotal,
      allResults[name]["details"][ruleKeyList[[ri]]],
      "with fallback"];
    If[TrueQ[showNoFallback],
      nfKey = Switch[ruleKeyList[[ri]], "R1", "RS", "R3", "RZ", _, Nothing];
      If[StringQ[nfKey],
        appendNetRowsArch[vers, faces, nf, nTotal,
          nfResults[name]["details"][nfKey],
          "no fallback"]]],
    {ri, Length[ruleKeyList]}
  ];

  (* SnubCube は鏡像（右手系）に異なる成功軌道がある（sq,snub 軌道）
     → 鏡像の結果も表示して全軌道タイプをカバーする *)
  If[name === "SnubCube",
    versMirror = vers * ConstantArray[{-1., 1., 1.}, Length[vers]];
    AppendTo[output,
      Style["  \[LongDash]\[LongDash] Mirror (right-handed Snub Cube) \[LongDash]\[LongDash]",
            13, GrayLevel[0.45]]];
    Do[
      ruleStr = Switch[ruleKeyList[[ri]], "R1", "maxphi", "R2", "loxo", "R3", "maxz", _, $Failed];
      If[ruleStr =!= $Failed,
        AppendTo[output,
          Style["  " <> ruleLabelList[[ri]], 13, ruleColorList[[ri]]]];
        mirrorRes = p3lRunAll[versMirror, faces, ruleStr];
        appendNetRowsArch[versMirror, faces, nf, nTotal, mirrorRes, "mirror"]
      ],
      {ri, Length[ruleKeyList]}
    ]
  ],
  {name, solidNames}
];

(* ----------------------------------------------------------------
   カラーバー
   ---------------------------------------------------------------- *)
AppendTo[output, Style[StringRepeat["\[LongDash]", 55], Gray]];
AppendTo[output,
  Style["\:30ab\:30e9\:30fc\:30d0\:30fc: \:9752=\:5148\:982d\:306e\:9762  \:7dd1=\:4e2d\:9593  \:8d64=\:672b\:5c3e\:306e\:9762  \
(\:5931\:6557\:6642\:306f\:8d64 \[Times] + \:914d\:7f6e\:6570/\:5168\:9762\:6570\:30e9\:30d9\:30eb)",
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
