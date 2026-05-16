(* ================================================================
   gen_120cell_faceup_nets.m
   120胞体 C1=3 の全隣接 C2（12通り）について
   A: cell-centroid-up  vs  B: 3D-face-centroid-up
   の3D展開図を並べて表示・PDF出力する．
   peel1Pair には peeling4Df4.m v4 の k>=3 Det=0 フォールバック修正を適用済み．
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName], NotebookDirectory[]];
Get[scriptDir <> "unfold3DExport.m"];

dataDir = scriptDir <> "_4DData/";
raw   = Get[dataDir <> "f120.m"];
vers  = N[raw[[1]]]; faces = raw[[4]]; cells = raw[[6]];
nc    = Length[cells];

fixedC1  = 3;
allC2    = {1, 2, 4, 10, 15, 17, 18, 19, 51, 97, 100, 115};

cellVerts[ci_] := Union @@ faces[[ cells[[ci]] ]];
wHat = {0., 0., 0., 1.};

(* ----------------------------------------------------------------
   face-up 回転行列
   ---------------------------------------------------------------- *)
faceUpRot[dir_] := Module[{d = Normalize[dir - Mean[vers]], ang},
  Which[
    Abs[d . wHat - 1] < 0.0001,  IdentityMatrix[4],
    Abs[d . wHat + 1] < 0.0001, -IdentityMatrix[4],
    True,
      ang = ArcCos[Clip[d . wHat, {-1, 1}]];
      RotationMatrix[ang, {d, wHat}]
  ]
];

(* ----------------------------------------------------------------
   RZ ペーリング（rot 方向を引数で指定）
   peeling4Df4.m v4 準拠：k>=3 Det=0 のとき xy クロス積フォールバック
   ---------------------------------------------------------------- *)
peel1Pair[c1_, c2_, rotDir_] :=
  Module[{eps = 10^-10, detEps = 10^-10, localVers, R, cC, c1Vec, c2Vec,
          adjList, order, last, newCell, cans, rCans,
          wvals, maxw, gvals, xyvals, tied, geoScore, xyScore},
    localVers = # - Mean[vers] & /@ vers;
    R         = faceUpRot[rotDir];
    localVers = R . # & /@ localVers;
    cC        = Mean[localVers[[ cellVerts[#] ]]] & /@ Range[nc];
    c1Vec = cC[[c1]]; c2Vec = cC[[c2]];
    adjList = Table[
      Select[Range[nc], j |-> j != i &&
        Length[Intersection[cells[[i]], cells[[j]]]] >= 1],
      {i, nc}];
    adjList = Complement[#, {c1}] & /@ adjList;
    order = {c1, c2};
    adjList = Complement[#, {c2}] & /@ adjList;
    last = c2;
    geoScore[k_] := If[Length[order] == 2,
      c2Vec[[1]] * cC[[k, 2]] - c2Vec[[2]] * cC[[k, 1]],
      Quiet[Det[{c1Vec, c2Vec, cC[[last]], cC[[k]]}]]];
    xyScore[k_] := c2Vec[[1]] * cC[[k, 2]] - c2Vec[[2]] * cC[[k, 1]];
    While[Length[Flatten[adjList]] != 0 && Length[adjList[[last]]] != 0,
      cans = If[Length[order] == 2, adjList[[last]],
        Select[adjList[[last]], j |->
          Quiet[Det[{c1Vec, c2Vec, cC[[last]], cC[[j]]}]] >= -eps]];
      newCell = If[Length[cans] >= 1,
        wvals = cC[[cans, 4]]; maxw = Max[wvals];
        tied  = Select[cans, cC[[#, 4]] >= maxw - 10^-10 &];
        If[Length[tied] == 1, tied[[1]],
          gvals = geoScore[#] & /@ tied;
          If[Length[order] >= 3 && Max[Abs[gvals]] < detEps,
            (* k>=3 Det=0 フォールバック: xy クロス積 *)
            xyvals = xyScore[#] & /@ tied;
            tied[[ First[Ordering[-xyvals, 1]] ]],
            tied[[ First[Ordering[-gvals, 1]] ]]]],
        rCans = adjList[[last]];
        wvals = cC[[rCans, 4]]; maxw = Min[wvals];
        tied  = Select[rCans, cC[[#, 4]] <= maxw + 10^-10 &];
        If[Length[tied] == 1, tied[[1]],
          gvals = Abs[geoScore[#]] & /@ tied;
          If[Length[order] >= 3 && Max[gvals] < detEps,
            xyvals = Abs[xyScore[#]] & /@ tied;
            tied[[ First[Ordering[-xyvals, 1]] ]],
            tied[[ First[Ordering[-gvals, 1]] ]]]]];
      last = newCell;
      AppendTo[order, newCell];
      adjList = Complement[#, {newCell}] & /@ adjList;
    ];
    {order, Length[order] == nc}
  ];

(* ----------------------------------------------------------------
   グレースケールグラデーション色: 0.5 (最初) → 0.9 (最後)
   ---------------------------------------------------------------- *)
apGray[k_, nTotal_] := GrayLevel[0.5 + 0.4 * (k - 1) / Max[nTotal - 1, 1]];

(* ----------------------------------------------------------------
   3D 展開図 Graphics3D を生成
   ---------------------------------------------------------------- *)
makeNet3D[ord_, embs_, sz_: 300] :=
  Module[{nC = Length[ord], prims},
    prims = Flatten[Table[
      Module[{ci = ord[[idx]], emb = embs[[idx]], col = apGray[idx, nC]},
        Map[Function[fi,
              {EdgeForm[{Thin, Black}],
               FaceForm[Directive[col, Opacity[0.85]]],
               Polygon[emb /@ faces[[fi]]]}],
            cells[[ci]]]],
      {idx, nC}], 1];
    Graphics3D[prims, Boxed -> False, ImageSize -> sz,
      Lighting -> "Neutral",
      ViewPoint -> {2.5, -2.5, 2.0}, SphericalRegion -> True]
  ];

(* ----------------------------------------------------------------
   precompute cellLocal3Ds
   ---------------------------------------------------------------- *)
Print["Pre-computing cellLocal3Ds ..."];
t0 = AbsoluteTime[];
cellLocal3Ds = Table[unfoldCellLocal3D[ci, vers, cells, faces], {ci, nc}];
Print["  Done: ", Round[AbsoluteTime[] - t0, 0.1], " s"];

c1CentroidDir = Mean[vers[[ cellVerts[fixedC1] ]]];

(* ----------------------------------------------------------------
   全 C2 について A・B を計算して Graphics3D を収集
   ---------------------------------------------------------------- *)
rows = {};
Do[
  c2 = allC2[[idx]];
  Print[Style["━━━ C2 = " <> ToString[c2] <> " ━━━", Bold, 14, GrayLevel[0.2]]];

  (* A: cell-centroid-up *)
  {ordA, okA} = peel1Pair[fixedC1, c2, c1CentroidDir];
  If[okA,
    embsA = unfoldTo3DFast[cells, faces, ordA, cellLocal3Ds];
    ovsA  = Length[checkOverlaps[ordA, cells, faces, embsA]],
    embsA = {}; ovsA = "—"
  ];

  (* B: 3D-face-centroid-up *)
  sharedFI = Intersection[cells[[fixedC1]], cells[[c2]]][[1]];
  sharedFC = Mean[vers[[ faces[[sharedFI]] ]]];
  {ordB, okB} = peel1Pair[fixedC1, c2, sharedFC];
  If[okB,
    embsB = unfoldTo3DFast[cells, faces, ordB, cellLocal3Ds];
    ovsB  = Length[checkOverlaps[ordB, cells, faces, embsB]],
    embsB = {}; ovsB = "—"
  ];

  Print["  A:", If[okA,"OK","FAIL"], " ov=", ovsA,
        "  B:", If[okB,"OK","FAIL"], " ov=", ovsB];

  (* ノートブック inline 表示 *)
  If[okA,
    Print[Style["  A: cell-centroid-up  (ov=" <> ToString[ovsA] <> ")", 11, GrayLevel[0.3]]];
    Print[makeNet3D[ordA, embsA, 350]]
  ];
  If[okB,
    Print[Style["  B: 3D-face-centroid-up  (ov=" <> ToString[ovsB] <> ")", 11, RGBColor[0.15, 0.35, 0.65]]];
    Print[makeNet3D[ordB, embsB, 350]]
  ];

  (* PDF 用 Grid 行 *)
  labelA = If[okA, "cell-up  (ov=" <> ToString[ovsA] <> ")", "cell-up  FAIL"];
  labelB = If[okB, "face-up  (ov=" <> ToString[ovsB] <> ")", "face-up  FAIL"];
  gA = If[okA, Rasterize[makeNet3D[ordA, embsA, 240], ImageResolution -> 150],
          Graphics[{Text[Style["FAIL", 18, Gray]]}, ImageSize -> {240, 240}]];
  gB = If[okB, Rasterize[makeNet3D[ordB, embsB, 240], ImageResolution -> 150],
          Graphics[{Text[Style["FAIL", 18, Gray]]}, ImageSize -> {240, 240}]];

  AppendTo[rows, {
    Style["C2=" <> ToString[c2], 11, Bold],
    Style[labelA, 10, GrayLevel[0.3]],
    Style[labelB, 10, GrayLevel[0.3]]
  }];
  AppendTo[rows, {"", gA, gB}],
  {idx, Length[allC2]}
];

(* ----------------------------------------------------------------
   出力
   ---------------------------------------------------------------- *)
header = {
  Style["", 1],
  Style["A: cell-centroid-up", Bold, 13, GrayLevel[0.2]],
  Style["B: 3D-face-centroid-up", Bold, 13, RGBColor[0.15, 0.35, 0.65]]
};

fig = Grid[
  Prepend[rows, header],
  Spacings -> {1, 0.4},
  Alignment -> {Center, Center},
  Background -> White
];

outFile = scriptDir <> "120cell_faceup_nets.pdf";
Export[outFile, Rasterize[fig, ImageResolution -> 200, Background -> White]];
Print["Exported: ", outFile];

(* ----------------------------------------------------------------
   3D-face-centroid-up 成功7ケースのみ Graphics3D で表示
   ---------------------------------------------------------------- *)
Print[Style["━━━ 3D-face-centroid-up: 成功ケース（7種）━━━", Bold, 14, RGBColor[0.15, 0.35, 0.65]]];
faceUpSuccC2 = {1, 2, 18, 19, 51, 97, 100};
Do[
  c2 = faceUpSuccC2[[i]];
  sharedFI = Intersection[cells[[fixedC1]], cells[[c2]]][[1]];
  sharedFC = Mean[vers[[ faces[[sharedFI]] ]]];
  {ordB, okB} = peel1Pair[fixedC1, c2, sharedFC];
  embsB = unfoldTo3DFast[cells, faces, ordB, cellLocal3Ds];
  ovsB  = Length[checkOverlaps[ordB, cells, faces, embsB]];
  Print[Style["  C2=" <> ToString[c2] <> "  (ov=" <> ToString[ovsB] <> ")",
              12, RGBColor[0.15, 0.35, 0.65]]];
  Print[makeNet3D[ordB, embsB, 400]],
  {i, Length[faceUpSuccC2]}];

(* ----------------------------------------------------------------
   STL 出力：face-up 成功7ケース × cell-up / face-up 両方
   ---------------------------------------------------------------- *)
stlDir = scriptDir <> "STLs_120cell_faceup/";
If[!DirectoryQ[stlDir], CreateDirectory[stlDir]];
Print[Style["━━━ STL 出力 ━━━", Bold, 14, GrayLevel[0.2]]];

exportEmbSTL[ord_, embs_, filename_] :=
  Module[{tris, gr},
    tris = unfoldTriangles[ord, cells, faces, embs];
    gr = Graphics3D[{EdgeForm[], Polygon /@ tris}];
    Export[filename, gr, "STL"];
    Print["  ", filename]];

Do[
  c2 = faceUpSuccC2[[i]];
  base = "120cell_C1-3_C2-" <> ToString[c2];

  (* A: cell-centroid-up *)
  {ordA, okA} = peel1Pair[fixedC1, c2, c1CentroidDir];
  embsA = unfoldTo3DFast[cells, faces, ordA, cellLocal3Ds];
  exportEmbSTL[ordA, embsA, stlDir <> base <> "_cellup.stl"];

  (* B: 3D-face-centroid-up *)
  sharedFI = Intersection[cells[[fixedC1]], cells[[c2]]][[1]];
  sharedFC = Mean[vers[[ faces[[sharedFI]] ]]];
  {ordB, okB} = peel1Pair[fixedC1, c2, sharedFC];
  embsB = unfoldTo3DFast[cells, faces, ordB, cellLocal3Ds];
  exportEmbSTL[ordB, embsB, stlDir <> base <> "_faceup.stl"],
  {i, Length[faceUpSuccC2]}];

Print["STL output complete: ", stlDir];
