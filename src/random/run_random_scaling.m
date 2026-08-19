(* run_random_scaling.m
   Apple-peel 成功率の N 依存性を調べる.
   N 個の母点から作った Delaunay / Voronoi 双対多面体に
   RS / RZ (fallback あり) を適用し, 成功率を N の関数としてプロット.
*)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
Get[baseDir <> "src/peel3d/peeling3DLoxo.m"];

sphericalRandom[] := Normalize[RandomVariate[NormalDistribution[], 3]];
dualToData[poly_] := {N[poly[[1]]], poly[[2]]};
successRate[results_] :=
  N[Length[Select[results, #[[4]] === True &]] / Length[results]];

(* ================================================================
   スケーリング計算
   N ごとにトライアル数を変える（大きい N は時間がかかるため）
   ================================================================ *)

nVals      = {20, 30, 50, 100, 200};
nTrialsVec = { 30, 20, 10,   5,   3};

Print["N依存性スケーリング計算開始"];
Print["N=", nVals, "  トライアル数=", nTrialsVec];

(* 結果格納: {N, DelRS, DelRZ, VorRS, VorRZ} の平均成功率 *)
scalingData = Table[
  Module[{nT = nTrialsVec[[ni]], rates},
    Print["\n--- N=", nVals[[ni]], "  (", nT, " trials) ---"];
    rates = Table[
      Module[{ptsT, mT, pT, vD, fD, vV, fV},
        ptsT = Table[sphericalRandom[], {nVals[[ni]]}];
        mT   = ConvexHullMesh[ptsT];
        pT   = DualPolyhedron[mT];
        vD   = MeshCoordinates[mT];
        fD   = MeshCells[mT, 2][[All, 1]];
        {vV, fV} = dualToData[pT];
        {successRate[p3lRunAll[vD, fD, "maxphi"]],  (* Del RS *)
         successRate[p3lRunAll[vD, fD, "maxz"]],    (* Del RZ *)
         successRate[p3lRunAll[vV, fV, "maxphi"]],  (* Vor RS *)
         successRate[p3lRunAll[vV, fV, "maxz"]]}    (* Vor RZ *)
      ],
      {nT}];
    With[{m = N[Mean[rates]], s = N[StandardDeviation[rates]]},
      Print["  Del RS=", Round[100 m[[1]], 0.1], "%",
            "  Del RZ=", Round[100 m[[2]], 0.1], "%",
            "  Vor RS=", Round[100 m[[3]], 0.1], "%",
            "  Vor RZ=", Round[100 m[[4]], 0.1], "%"];
      Print["  (std: Del RS=", Round[100 s[[1]], 0.1],
            "  Del RZ=", Round[100 s[[2]], 0.1],
            "  Vor RS=", Round[100 s[[3]], 0.1],
            "  Vor RZ=", Round[100 s[[4]], 0.1], ")"];
      Prepend[m, nVals[[ni]]]]],
  {ni, Length[nVals]}];

(* ================================================================
   サマリ表
   ================================================================ *)

Print["\n", StringRepeat["=", 60]];
Print["N依存性サマリ (平均成功率 %)"];
Print[StringRepeat["=", 60]];
Print[TableForm[
  Prepend[
    {ToString[#[[1]]], Round[100 #[[2]], 0.1],
     Round[100 #[[3]], 0.1], Round[100 #[[4]], 0.1],
     Round[100 #[[5]], 0.1]} & /@ scalingData,
    {"N", "Del RS", "Del RZ", "Vor RS", "Vor RZ"}],
  TableSpacing -> {1, 2}]];

(* 面次数の平均の理論値: (6N-12)/N = 6 - 12/N *)
Print["\nVoronoi 平均面次数の理論値 (6 - 12/N):"];
Do[Print["  N=", n, ": ", N[6 - 12/n, 3]], {n, nVals}];

(* ================================================================
   プロット
   ================================================================ *)

nPlot   = scalingData[[All, 1]];
delRZ   = 100. scalingData[[All, 3]];
vorRZ   = 100. scalingData[[All, 5]];
delRS   = 100. scalingData[[All, 2]];
vorRS   = 100. scalingData[[All, 4]];

plt = ListLinePlot[
  {Transpose[{nPlot, vorRZ}],
   Transpose[{nPlot, vorRS}],
   Transpose[{nPlot, delRZ}],
   Transpose[{nPlot, delRS}]},
  PlotMarkers -> {Automatic, 10},
  PlotStyle  -> {
    {Blue,  Thickness[0.003]},
    {Cyan,  Thickness[0.003]},
    {Red,   Thickness[0.003]},
    {Pink,  Thickness[0.003]}},
  PlotLegends -> Placed[
    {"Voronoi RZ", "Voronoi RS", "Delaunay RZ", "Delaunay RS"},
    {0.75, 0.6}],
  AxesLabel  -> {"N (母点数)", "成功率 (%)"},
  PlotLabel  -> Style["Apple-Peel 成功率 vs N\n(Delaunay vs Voronoi 双対)", 13, Bold],
  PlotRange  -> {{15, 210}, {-2, 102}},
  GridLines  -> Automatic,
  GridLinesStyle -> LightGray,
  ImageSize  -> 500];

Print[plt];
Export[baseDir <> "random_scaling_plot.png", plt, ImageSize -> 700];
Print["プロット保存: random_scaling_plot.png"];

(* 参照: 正多面体・アルキメデス多面体の成功率 *)
Print["\n--- 参照: 既知の結果 ---"];
Print["正多面体 (全規則, fallback あり): 全5種 100%"];
Print["アルキメデス RZ: 三角形なし6種は高成功率, 三角形多数は低成功率"];
Print["4D: Voronoi 類似(8/120-cell)=Perfect, Delaunay 類似(16/600-cell)=Possible/Impossible"];
