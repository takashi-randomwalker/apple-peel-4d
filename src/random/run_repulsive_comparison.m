(* run_repulsive_comparison.m
   Thomson 反発エネルギー最小化（Coulomb 勾配降下法）で
   点分布をランダムから均一へ変化させ，
   Voronoi 双対多面体の Apple-Peel RZ 成功率の変化を調べる.

   パラメータ化:
     relaxLevel = 0  : 純粋ランダム（Poisson 点過程）
     relaxLevel = ∞  : Thomson 平衡点（Coulomb エネルギー最小）
   勾配降下ステップ数で補間する.
*)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
Get[baseDir <> "src/peel3d/peeling3DLoxo.m"];

sphericalRandom[] := Normalize[RandomVariate[NormalDistribution[], 3]];
dualToData[poly_] := {N[poly[[1]]], poly[[2]]};
successRate[results_] :=
  N[Length[Select[results, #[[4]] === True &]] / Length[results]];

(* ================================================================
   Thomson 勾配降下法
   Coulomb エネルギー E = Σ_{i<j} 1/||p_i - p_j|| を最小化.
   各ステップ: 勾配方向に動かし → 球面に投影（正規化）.
   ================================================================ *)

(* 各点 i に働く Coulomb 斥力の勾配（球面接平面成分） *)
coulombGradient[pts_] :=
  Module[{n = Length[pts]},
    Table[
      Module[{gi},
        gi = Sum[
          With[{d = pts[[i]] - pts[[j]]},
            If[i === j, {0.,0.,0.}, d / Norm[d]^3]],
          {j, n}];
        (* 球面接平面への射影: g - (g·p)p *)
        gi - (gi . pts[[i]]) pts[[i]]],
      {i, n}]];

(* nSteps ステップの勾配降下で緩和した点を返す *)
relaxPoints[pts0_, nSteps_Integer, lr_Real: 0.005] :=
  Nest[
    Function[pts,
      Normalize /@ (pts - lr * coulombGradient[pts])],
    N[pts0],
    nSteps];

(* Coulomb エネルギー（均一性の指標）*)
coulombEnergy[pts_] :=
  N[Sum[1/Norm[pts[[i]] - pts[[j]]], {i, Length[pts]}, {j, i+1, Length[pts]}]];

(* Voronoi 面次数の分散（均一性の指標）: 0 に近いほど均一 *)
faceDegreeSd[facesVor_] := N[StandardDeviation[Length /@ facesVor]];

(* ================================================================
   実験 1: 固定 N=30, 緩和ステップ数を変化
   ================================================================ *)

Print["\n", StringRepeat["=", 60]];
Print["実験1: 緩和レベルと成功率 (N=30, 10 trials/level)"];
Print[StringRepeat["=", 60]];

nPts   = 30;
nTrial = 10;
stepLevels = {0, 10, 50, 200, 500, 2000};

relaxResults = Table[
  Module[{rates, energies, fsdList},
    {rates, energies, fsdList} = Transpose[Table[
      Module[{pts0, pts, mT, pT, vV, fV, rVRZ},
        pts0 = Table[sphericalRandom[], {nPts}];
        pts  = relaxPoints[pts0, nSteps];
        mT   = ConvexHullMesh[pts];
        pT   = DualPolyhedron[mT];
        {vV, fV} = dualToData[pT];
        rVRZ = p3lRunAll[vV, fV, "maxz"];
        {successRate[rVRZ],
         coulombEnergy[pts],
         faceDegreeSd[fV]}],
      {nTrial}]];
    With[{mr = Mean[rates], me = Mean[energies], mf = Mean[fsdList]},
      Print["  steps=", PaddedForm[nSteps, 5],
            "  VorRZ=", PaddedForm[Round[100. mr, 0.1], {5,1}], "%",
            "  E_C=", PaddedForm[Round[me, 0.01], {7,2}],
            "  SD(degree)=", PaddedForm[Round[mf, 0.01], {4,2}]];
      {nSteps, mr, me, mf}]],
  {nSteps, stepLevels}];

(* ================================================================
   実験 2: N=20, 30, 50 × ランダム vs Thomson 平衡
   ================================================================ *)

Print["\n", StringRepeat["=", 60]];
Print["実験2: ランダム vs Thomson 平衡  (各 10 trials)"];
Print[StringRepeat["=", 60]];
Print["  各 N について: steps=0 (ランダム) vs steps=2000 (Thomson 近似)"];

nTrialComp = 10;
stepsEquil = 2000;

Do[
  Module[{rRand, rRelax, eRand, eRelax},
    {rRand, eRand} = Transpose[Table[
      Module[{pts, mT, pT, vV, fV, r},
        pts = Table[sphericalRandom[], {nN}];
        mT  = ConvexHullMesh[pts];
        pT  = DualPolyhedron[mT];
        {vV, fV} = dualToData[pT];
        r = p3lRunAll[vV, fV, "maxz"];
        {successRate[r], coulombEnergy[pts]}],
      {nTrialComp}]];
    {rRelax, eRelax} = Transpose[Table[
      Module[{pts0, pts, mT, pT, vV, fV, r},
        pts0 = Table[sphericalRandom[], {nN}];
        pts  = relaxPoints[pts0, stepsEquil];
        mT   = ConvexHullMesh[pts];
        pT   = DualPolyhedron[mT];
        {vV, fV} = dualToData[pT];
        r = p3lRunAll[vV, fV, "maxz"];
        {successRate[r], coulombEnergy[pts]}],
      {nTrialComp}]];
    Print["\n  N=", nN];
    Print["    random :  VorRZ=", Round[100. Mean[rRand], 0.1],
          "%  E_C=", Round[Mean[eRand], 1]];
    Print["    Thomson:  VorRZ=", Round[100. Mean[rRelax], 0.1],
          "%  E_C=", Round[Mean[eRelax], 1]]],
  {nN, {20, 30, 50}}];

(* ================================================================
   プロット: 緩和ステップ vs 成功率・Coulomb エネルギー
   ================================================================ *)

steps = relaxResults[[All, 1]];
rates = 100. relaxResults[[All, 2]];
eners = relaxResults[[All, 3]];

pltRate = ListLinePlot[
  Transpose[{steps, rates}],
  PlotMarkers -> {Automatic, 10},
  PlotStyle   -> {Blue, Thickness[0.003]},
  AxesLabel   -> {"緩和ステップ数", "VorRZ 成功率 (%)"},
  PlotLabel   -> Style["Thomson 緩和 vs Apple-Peel 成功率\n(N=30, RZ, 10 trials/level)", 12, Bold],
  PlotRange   -> {{-50, 2100}, {-3, 103}},
  GridLines   -> Automatic, GridLinesStyle -> LightGray,
  ImageSize   -> 420];

pltEnergy = ListLinePlot[
  Transpose[{steps, eners}],
  PlotMarkers -> {Automatic, 10},
  PlotStyle   -> {Red, Thickness[0.003]},
  AxesLabel   -> {"緩和ステップ数", "Coulomb エネルギー"},
  PlotLabel   -> Style["Thomson 緩和 vs Coulomb エネルギー\n(N=30, 10 trials/level)", 12, Bold],
  PlotRange   -> {{-50, 2100}, Automatic},
  GridLines   -> Automatic, GridLinesStyle -> LightGray,
  ImageSize   -> 420];

combined = GraphicsRow[{pltRate, pltEnergy}, Spacings -> 20];
Print[combined];
Export[baseDir <> "repulsive_comparison.png", combined, ImageSize -> 1000];
Print["プロット保存: repulsive_comparison.png"];
