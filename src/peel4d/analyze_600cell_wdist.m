(* ================================================================
   analyze_600cell_wdist.m
   ================================================================
   600胞体のセル重心 w 分布（cell-centroid-up 後）を計算し，
   Kaino (2019) の vertex-first 層構造と比較する．

   Kaino の 13 層: [20,20,30,60,60,80,60,80,60,60,30,20,20]
   累積境界: 20, 40, 70, 130, 190, 270, 330, 410, 470, 530, 560, 580, 600
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];
dataDir = scriptDir <> "_4DData/";

Print["Loading 600-cell data..."];
raw = Get[dataDir <> "f600.m"];
vers  = N[raw[[1]]];
faces = raw[[4]];
cells = raw[[6]];
nc    = Length[cells];
Print["Cells: ", nc, "  Vertices: ", Length[vers]];

(* セル頂点リスト *)
cellVerts[i_] := Union @@ faces[[ cells[[i]] ]];

(* C1 = top (代表として 1 を使用) *)
top = 1;

(* Step 1: 重心を原点へ *)
localVers = (# - Mean[vers]) & /@ vers;

(* Step 2: C1 の重心 *)
cc1 = Mean[ localVers[[ cellVerts[top] ]] ];

(* Step 3: cell-centroid-up (C1 重心を +w 方向へ) *)
wHat = {0, 0, 0, 1};
localVers = Which[
  Abs[(cc1 / Norm[cc1]) . wHat - 1] < 0.0001,
    localVers,
  Abs[(cc1 / Norm[cc1]) . wHat + 1] < 0.0001,
    -localVers,
  True,
    angle4D = ArcCos[Clip[cc1 . wHat / Norm[cc1], {-1, 1}]];
    RotationMatrix[angle4D, {cc1 / Norm[cc1], wHat}] . # & /@ localVers
];

(* Step 4: 全セル重心を計算 *)
cCenters = Mean[ localVers[[ cellVerts[#] ]] ] & /@ Range[nc];
wCoords  = cCenters[[All, 4]];

Print["\n--- 600-cell cell centroid w-distribution (top=", top, ") ---"];
Print["w range: [", Min[wCoords] // N // Round[#, 0.001]&, ", ",
                     Max[wCoords] // N // Round[#, 0.001]&, "]"];

(* w 座標でクラスタリング（近い値をまとめる） *)
wSorted = Sort[wCoords];
eps = 0.001;

(* 連続値をグループ化 *)
groups = {};
current = {wSorted[[1]]};
Do[
  If[Abs[wSorted[[i]] - Last[current]] < eps,
    AppendTo[current, wSorted[[i]]],
    AppendTo[groups, current];
    current = {wSorted[[i]]}
  ],
  {i, 2, Length[wSorted]}
];
AppendTo[groups, current];

Print["\n--- Layer structure (", Length[groups], " bands) ---"];
Print[StringPadLeft["Layer", 6], "  ",
      StringPadLeft["Count", 6], "  ",
      StringPadLeft["w-value", 10], "  ",
      StringPadLeft["Cumulative", 11]];
Print[StringRepeat["-", 40]];

cumul = 0;
layerCounts = {};
layerWVals  = {};
Do[
  cnt = Length[groups[[i]]];
  wval = Mean[groups[[i]]];
  cumul += cnt;
  AppendTo[layerCounts, cnt];
  AppendTo[layerWVals, wval];
  Print[StringPadLeft[ToString[i], 6], "  ",
        StringPadLeft[ToString[cnt], 6], "  ",
        StringPadLeft[ToString[Round[wval, 0.0001]], 10], "  ",
        StringPadLeft[ToString[cumul], 11]],
  {i, Length[groups]}
];

Print["\n--- Kaino (2019) vertex-first comparison ---"];
kainoLayers = {20, 20, 30, 60, 60, 80, 60, 80, 60, 60, 30, 20, 20};
kainoCumul  = Accumulate[kainoLayers];
Print["Kaino layers: ", kainoLayers];
Print["Kaino cumul:  ", kainoCumul];
Print["Our   layers: ", layerCounts];
Print["Our   cumul:  ", Accumulate[layerCounts]];

(* 停止ステップとの対応 *)
Print["\n--- Stopping step analysis ---"];
stoppingSteps = {146, 150, 276, 279, 284};
Print["Observed stopping steps: ", stoppingSteps];
Print[""];
Print["At each stopping step, the algorithm has placed this many cells:"];
Do[
  s = stoppingSteps[[i]];
  (* 停止ステップ s: order に s 個セルが入っている = step s で詰まる *)
  (* 実際には step s で stuck になる = order length は s *)
  cumOur  = Accumulate[layerCounts];
  (* どの層境界に最も近いか *)
  diffs   = Abs[cumOur - s];
  nearest = First[Ordering[diffs, 1]];
  Print["  step ", s, ": nearest our-layer boundary = ", cumOur[[nearest]],
        " (layer ", nearest, "/", Length[groups], ")",
        "  diff = ", diffs[[nearest]]],
  {i, Length[stoppingSteps]}
];

Print["\n--- Kaino cumulative boundaries near stopping steps ---"];
Do[
  s = stoppingSteps[[i]];
  diffs   = Abs[kainoCumul - s];
  nearest = First[Ordering[diffs, 1]];
  Print["  step ", s, ": nearest Kaino boundary = ", kainoCumul[[nearest]],
        " (after layer ", nearest, ")",
        "  diff = ", diffs[[nearest]]],
  {i, Length[stoppingSteps]}
];

(* C2 の w 座標確認（k=2 タイブレーク関連） *)
Print["\n--- C1=", top, " neighbors (C2 candidates) ---"];
neighborsC1 = Select[
  Range[nc],
  Function[j, j != top && Length[Intersection[cells[[top]], cells[[j]]]] >= 1]
];
Print["Number of C2 candidates: ", Length[neighborsC1]];
c2wVals = cCenters[[neighborsC1, 4]];
Print["C2 w-values: ", Round[c2wVals, 0.0001]];
Print["C2 w-spread: ", Max[c2wVals] - Min[c2wVals]];
Print["C1 w-value:  ", Round[cCenters[[top, 4]], 0.0001]];

Print["\nDone."];
