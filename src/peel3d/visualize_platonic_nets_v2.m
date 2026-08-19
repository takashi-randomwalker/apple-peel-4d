(* ================================================================
   visualize_platonic_nets_v2.m
   Apple-Peel 展開図の可視化 (現行版)

   変更点 (v2):
     - R2 (loxodrome) を廃止; RS (Spiral rule) と RZ (Zonal rule) のみ表示
     - 規則名を RS / RZ に統一 (R1 / R3 → RS / RZ)
     - dataPlatonic_strict.mx (eps < 0, thresh = +10^-10) の表示を追加
       showStrict = True  : 標準 (eps > 0) + 厳密 (eps < 0) の両方を表示
       showStrict = False : 標準 (eps > 0) のみ表示

   データソース:
     dataPlatonic.mx       ... RS / RZ (標準, eps > 0, thresh = -10^-10)
     dataPlatonic_strict.mx ... RS / RZ (厳密, eps < 0, thresh = +10^-10)

   展開アルゴリズム:
     1. F1 の外向き法線が +z 軸になるよう全頂点を回転 (F1 が xy 平面に乗る)
     2. k=2..n: F(k-1) と Fk の共有辺を軸に Fk を xy 平面に倒す
     3. c1→c2 が +x 方向になるよう最終回転

   使い方:
     Get["/path/to/visualize_platonic_nets_v2.m"]
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
If[scriptDir === "", scriptDir = NotebookDirectory[]];
Get[scriptDir <> "peeling3DLoxo.m"];
Get[scriptDir <> "dataPlatonic.mx"];
Get[scriptDir <> "dataPlatonic_strict.mx"];

(* ----------------------------------------------------------------
   設定パラメータ
   ---------------------------------------------------------------- *)
If[!ValueQ[netSize],        netSize        = 120];       (* 各展開図の画像サイズ (px) *)
If[!ValueQ[maxShowOK],      maxShowOK      = Infinity];  (* 成功展開図の表示上限 *)
If[!ValueQ[maxShowFail],    maxShowFail    = Infinity];  (* 失敗展開図の表示上限 *)
If[!ValueQ[deduplicate],    deduplicate    = True];      (* True: 同じ order を1件に集約 *)
If[!ValueQ[showStrict],     showStrict     = True];      (* True: 厳密条件 (eps < 0) も表示 *)
If[!ValueQ[deduplicateGeo], deduplicateGeo = False];     (* True: 幾何学的重複除去（鏡像含む） *)

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
   幾何学的重複除去キー: 面重心のペアワイズ距離（ソート済み）
   回転・並進・スケール・反射不変 → 鏡像も同一キーになる
   ---------------------------------------------------------------- *)
netGeoKey[vers_, faces_, order_] :=
  Module[{polys2D, centroids, c1, shifted, maxD},
    polys2D   = p3lUnfoldNet[vers, faces, order];
    centroids = Mean /@ polys2D;
    c1        = centroids[[1]];
    shifted   = # - c1 & /@ centroids;
    maxD      = Max[Norm /@ shifted];
    If[maxD < 10^-8, Return[{}]];
    shifted = shifted / maxD;
    Round[Sort@Flatten@Table[
      Norm[shifted[[i]] - shifted[[j]]],
      {i, Length[shifted]}, {j, i + 1, Length[shifted]}
    ], 10^-3]
  ];

(* ----------------------------------------------------------------
   グラデーション色: 青(先頭) -> 緑(中間) -> 赤(末尾)
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
   失敗展開図の描画 (赤 x マーク・配置数/全数ラベル)
   ---------------------------------------------------------------- *)
netGraphicsFail[polys2D_, nPlaced_, nfTotal_, sz_: 140] :=
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
   1変種（rule × fallback）の展開図行を output に追記するヘルパー
   ---------------------------------------------------------------- *)
appendNetRows[vers_, faces_, nf_, nTotal_, allRes_, label_] :=
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

    AppendTo[output,
      Style["    [" <> label <> "]  \:6210\:529f: " <> ToString[nComp] <>
            "/" <> ToString[nTotal] <>
            "   \:5931\:6557: " <> ToString[nTotal - nComp] <> "/" <> ToString[nTotal] <>
            If[deduplicateGeo, "  (geo-uniq: " <> ToString[Length[uniqOK]] <> ")", ""],
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
   メインループ
   ================================================================ *)
solidNames  = {"Tetrahedron", "Cube", "Octahedron",
               "Dodecahedron", "Icosahedron"};

(* 標準 (eps > 0): dataPlatonic.mx では R1 / R3 キーを使用 *)
stdRuleKeys   = {"R1", "R3"};
stdRuleLabels = {"RS (Spiral rule, max \[Phi])", "RZ (Zonal rule, max z)"};

(* 厳密 (eps < 0): dataPlatonic_strict.mx では RS / RZ キーを使用 *)
strRuleKeys   = {"RS", "RZ"};
strRuleLabels = {"RS (Spiral rule, max \[Phi])", "RZ (Zonal rule, max z)"};

fbKeys   = {"withFallback", "noFallback"};
fbLabels = {"with fallback", "no fallback"};

output = {};

(* ----------------------------------------------------------------
   標準条件 (eps > 0, thresh = -10^-10)
   ---------------------------------------------------------------- *)
AppendTo[output,
  Style["\:6a19\:6e96\:6761\:4ef6  (eps > 0,  thresh = \[Minus]10^\[Minus]10)  det \[GreaterEqual] \[Minus]\[Epsilon]",
        Bold, 14, RGBColor[0.1, 0.3, 0.6]]];

Do[
  {vers, faces} = p3lExtractPolyhedron[name];
  nf     = platonicData[name]["nFaces"];
  nTotal = platonicData[name]["totalPairs"];

  AppendTo[output, Style[StringRepeat["\[LongDash]", 50], Gray]];
  AppendTo[output,
    Style[name <> "  (faces: " <> ToString[nf] <> ")", Bold, 16]];

  Do[
    AppendTo[output,
      Style["  " <> stdRuleLabels[[ri]], 14, RGBColor[0.2, 0.2, 0.6]]];
    Do[
      appendNetRows[vers, faces, nf, nTotal,
        platonicData[name][stdRuleKeys[[ri]]][fbKeys[[fi]]],
        fbLabels[[fi]]],
      {fi, 2}
    ],
    {ri, 2}
  ],
  {name, solidNames}
];

(* ----------------------------------------------------------------
   厳密条件 (eps < 0, thresh = +10^-10)   (showStrict = True の場合)
   ---------------------------------------------------------------- *)
If[showStrict,
  AppendTo[output, Style[StringRepeat["\[LongDash]", 50], Gray]];
  AppendTo[output,
    Style["\:53b3\:5bc6\:6761\:4ef6  (eps < 0,  thresh = +10^\[Minus]10)  det \[GreaterEqual] +\[Epsilon]",
          Bold, 14, RGBColor[0.5, 0.1, 0.1]]];

  Do[
    {vers, faces} = p3lExtractPolyhedron[name];
    nf     = platonicDataStrict[name]["nFaces"];
    nTotal = platonicDataStrict[name]["totalPairs"];

    AppendTo[output, Style[StringRepeat["\[LongDash]", 50], Gray]];
    AppendTo[output,
      Style[name <> "  (faces: " <> ToString[nf] <> ")", Bold, 16]];

    Do[
      AppendTo[output,
        Style["  " <> strRuleLabels[[ri]], 14, RGBColor[0.5, 0.1, 0.1]]];
      Do[
        appendNetRows[vers, faces, nf, nTotal,
          platonicDataStrict[name][strRuleKeys[[ri]]][fbKeys[[fi]]],
          fbLabels[[fi]]],
        {fi, 2}
      ],
      {ri, 2}
    ],
    {name, solidNames}
  ]
];

(* ----------------------------------------------------------------
   カラーバー
   ---------------------------------------------------------------- *)
AppendTo[output, Style[StringRepeat["\[LongDash]", 50], Gray]];
AppendTo[output,
  Style["\:30ab\:30e9\:30fc\:30d0\:30fc: \:9752=\:5148\:982d\:306e\:9762  \:7dd1=\:4e2d\:9593  \:8d64=\:672b\:5c3e\:306e\:9762" <>
        "  (\:5931\:6557\:6642\:306f\:8d64 \[Times] + \:914d\:7f6e\:6570/\:5168\:9762\:6570\:30e9\:30d9\:30eb)",
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
