(* run_random_polyhedra.m
   Apple-peel RS/RZ を Delaunay 三角分割 vs Voronoi 双対多面体で比較する.
   同一の母点 pts から ConvexHullMesh (全三角形, 3-正則) と
   DualPolyhedron (混合多角形, 平均次数≈6) を生成し, 成功率を比較する.

   実行: Get[baseDir <> "src/random/run_random_polyhedra.m"]
*)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
Get[baseDir <> "src/peel3d/peeling3DLoxo.m"];

(* ================================================================
   ユーティリティ
   ================================================================ *)

sphericalRandom[] := Normalize[RandomVariate[NormalDistribution[], 3]];

(* DualPolyhedron (Polyhedron 式) からデータを取り出す *)
dualToData[poly_] := {N[poly[[1]]], poly[[2]]};

(* 面次数の分布を集計 *)
degreeStats[faces_] :=
  Module[{degs = Length /@ faces},
    <|"min"  -> Min[degs],
      "max"  -> Max[degs],
      "mean" -> N[Mean[degs]],
      "dist" -> Sort[Tally[degs]]|>];

(* 結果から成功率・失敗時ステップ数中央値を計算
   results の各要素: {top, f2i, order, completeBool} *)
analyzeResults[results_] :=
  Module[{nTotal, nOK, nFail, steps},
    nTotal = Length[results];
    nOK    = Length[Select[results, #[[4]] === True &]];
    nFail  = nTotal - nOK;
    steps  = Length[#[[3]]] & /@ Select[results, #[[4]] === False &];
    <|"total"    -> nTotal,
      "ok"       -> nOK,
      "rate"     -> If[nTotal > 0, N[nOK / nTotal], 0],
      "failStepMedian" ->
        If[nFail > 0, N[Median[steps]], Missing["NoFailure"]]|>];

(* ================================================================
   1インスタンス詳細解析
   ================================================================ *)

SeedRandom[42];
nPts = 30;
pts  = Table[sphericalRandom[], {nPts}];

m = ConvexHullMesh[pts];
p = DualPolyhedron[m];

(* データ抽出 *)
versDel  = MeshCoordinates[m];
facesDel = MeshCells[m, 2][[All, 1]];
{versVor, facesVor} = dualToData[p];

nFDel = Length[facesDel];
nFVor = Length[facesVor];

Print["\n", StringRepeat["=", 55]];
Print["ランダム多面体 Apple-Peel 比較  (N=", nPts, " 母点)"];
Print[StringRepeat["=", 55]];

Print["\n--- Delaunay 三角分割 ---"];
Print["  面数: ", nFDel, "  (= 2N-4 = ", 2 nPts - 4, ")"];
Print["  全面が三角形, 面隣接次数 = 3"];
Print["  (top,f2) ペア数: ", Total[Length /@ p3lBuildAdj[facesDel]]];

Print["\n--- Voronoi 双対多面体 ---"];
Print["  面数: ", nFVor, "  (= N = ", nPts, ")"];
With[{ds = degreeStats[facesVor]},
  Print["  面次数: min=", ds["min"], "  max=", ds["max"],
        "  mean=", Round[ds["mean"], 0.01]];
  Print["  分布: ", ds["dist"]];
  Print["  Σ(6-k) = ",
        Total[(6 - Length[#]) & /@ facesVor], "  (Euler, 期待値 12)"]];
Print["  (top,f2) ペア数: ", Total[Length /@ p3lBuildAdj[facesVor]],
      "  (双対なので Delaunay と同数)"];

(* ================================================================
   Apple-Peel 実行 (RS / RZ, fallback あり)
   ================================================================ *)

Print["\n", StringRepeat["-", 55]];
Print["Apple-Peel 実行中 (fallback あり)..."];

{tDelRS, resDelRS} = AbsoluteTiming[p3lRunAll[versDel, facesDel, "maxphi"]];
{tDelRZ, resDelRZ} = AbsoluteTiming[p3lRunAll[versDel, facesDel, "maxz"]];
{tVorRS, resVorRS} = AbsoluteTiming[p3lRunAll[versVor, facesVor, "maxphi"]];
{tVorRZ, resVorRZ} = AbsoluteTiming[p3lRunAll[versVor, facesVor, "maxz"]];

aDelRS = analyzeResults[resDelRS];
aDelRZ = analyzeResults[resDelRZ];
aVorRS = analyzeResults[resVorRS];
aVorRZ = analyzeResults[resVorRZ];

(* ================================================================
   結果サマリ
   ================================================================ *)

Print["\n", StringRepeat["=", 55]];
Print["結果サマリ (N=", nPts, ", seed=42)"];
Print[StringRepeat["=", 55]];

fmtRate[a_] := ToString[a["ok"]] <> "/" <> ToString[a["total"]] <>
               "  (" <> ToString[Round[100. a["rate"], 0.1]] <> "%)";
fmtStep[a_] := If[MissingQ[a["failStepMedian"]],
                  "—", ToString[a["failStepMedian"]]];

Print[TableForm[{
  {"", "RS", "RZ"},
  {"Delaunay (" <> ToString[nFDel] <> "面, 3-正則)",
   fmtRate[aDelRS], fmtRate[aDelRZ]},
  {"Voronoi  (" <> ToString[nFVor] <> "面, 平均" <>
   ToString[Round[Mean[Length /@ facesVor], 0.1]] <> "-正則)",
   fmtRate[aVorRS], fmtRate[aVorRZ]}
}, TableSpacing -> {1, 2}]];

Print["失敗時ステップ中央値 (全面数: Del=",nFDel,", Vor=",nFVor,"):"];
Print["  Delaunay  RS=", fmtStep[aDelRS], "  RZ=", fmtStep[aDelRZ]];
Print["  Voronoi   RS=", fmtStep[aVorRS], "  RZ=", fmtStep[aVorRZ]];

(* ================================================================
   複数インスタンスで統計 (nTrial 回)
   ================================================================ *)

Print["\n", StringRepeat["=", 55]];
nTrial = 20;
Print["複数インスタンス統計 (N=", nPts, ", ", nTrial, " trial)"];
Print[StringRepeat["=", 55]];

trialRates = Table[
  Module[{ptsT, mT, pT, vD, fD, vV, fV, rDRS, rDRZ, rVRS, rVRZ},
    ptsT = Table[sphericalRandom[], {nPts}];
    mT   = ConvexHullMesh[ptsT];
    pT   = DualPolyhedron[mT];
    vD   = MeshCoordinates[mT];
    fD   = MeshCells[mT, 2][[All, 1]];
    {vV, fV} = dualToData[pT];
    rDRS = p3lRunAll[vD, fD, "maxphi"];
    rDRZ = p3lRunAll[vD, fD, "maxz"];
    rVRS = p3lRunAll[vV, fV, "maxphi"];
    rVRZ = p3lRunAll[vV, fV, "maxz"];
    {analyzeResults[rDRS]["rate"],
     analyzeResults[rDRZ]["rate"],
     analyzeResults[rVRS]["rate"],
     analyzeResults[rVRZ]["rate"]}],
  {nTrial}];

meanRates = N[Mean[trialRates]];

Print[TableForm[{
  {"", "RS", "RZ"},
  {"Delaunay (平均成功率)",
   ToString[Round[100. meanRates[[1]], 0.1]] <> "%",
   ToString[Round[100. meanRates[[2]], 0.1]] <> "%"},
  {"Voronoi  (平均成功率)",
   ToString[Round[100. meanRates[[3]], 0.1]] <> "%",
   ToString[Round[100. meanRates[[4]], 0.1]] <> "%"}
}, TableSpacing -> {1, 2}]];

Print["\n4D 正多胞体との対比:"];
Print["  Delaunay (3-正則) ↔ 16-cell/600-cell (Possible/Impossible)"];
Print["  Voronoi  (6-正則) ↔ 8-cell/120-cell  (Perfect)"];
