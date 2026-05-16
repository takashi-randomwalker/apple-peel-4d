(* ================================================================
   analyze_c2_xyz.m
   120胞体 C1=3 の隣接 C2 重心の xyz 構造を分析
   ================================================================ *)
scriptDir = If[ValueQ[$InputFileName], DirectoryName[$InputFileName], NotebookDirectory[]];
dataDir = scriptDir <> "_4DData/";
raw = Get[dataDir <> "f120.m"];
vers = N[raw[[1]]]; faces = raw[[4]]; cells = raw[[6]];
nc = Length[cells]; fixedC1 = 3;
allC2 = {1, 2, 4, 10, 15, 17, 18, 19, 51, 97, 100, 115};
cellVerts[ci_] := Union @@ faces[[cells[[ci]]]];
wHat = {0., 0., 0., 1.};

(* cell-centroid-up 回転 *)
localVers = # - Mean[vers] & /@ vers;
c1d = Normalize[Mean[vers[[cellVerts[fixedC1]]]] - Mean[vers]];
R = Which[
  Abs[c1d . wHat - 1] < 0.0001, IdentityMatrix[4],
  Abs[c1d . wHat + 1] < 0.0001, -IdentityMatrix[4],
  True, RotationMatrix[ArcCos[Clip[c1d . wHat, {-1, 1}]], {c1d, wHat}]];
rotVers = R . # & /@ localVers;
cC = Mean[rotVers[[cellVerts[#]]]] & /@ Range[nc];

(* C2 重心テーブル（全座標） *)
Print[Style["=== 12 C2 centroids in cell-centroid-up frame ===", Bold]];
tab = Table[
  pt = cC[[c2]];
  rxy  = Sqrt[pt[[1]]^2 + pt[[2]]^2];
  rxyz = Norm[pt[[{1,2,3}]]];
  th   = ArcTan[pt[[1]], pt[[2]]] * 180/Pi;
  {c2, Round[pt[[1]],0.001], Round[pt[[2]],0.001], Round[pt[[3]],0.001],
   Round[pt[[4]],4], Round[rxy,0.001], Round[rxyz,0.001], Round[th,0.1]},
  {c2, allC2}];
Print[TableForm[tab,
  TableHeadings -> {None, {"C2","x","y","z","w","r_xy","r_xyz","θ_xy(°)"}}]];

(* z 値の分布 *)
zVals = cC[[allC2, 3]];
Print["\nz values: ", Round[Sort[zVals, Greater], 0.001]];
Print["Unique z values: ", Round[Sort[Union[Round[zVals, 0.001]], Greater], 0.001]];
Print["r_xyz values: ", Round[Sort[Norm /@ cC[[allC2, {1,2,3}]], Greater], 0.001]];

(* ペアワイズ距離 *)
xyz = cC[[allC2, {1,2,3}]];
pdists = Flatten[Table[Norm[xyz[[i]] - xyz[[j]]], {i, 12}, {j, i+1, 12}]];
udists = Sort[Union[Round[pdists, 0.001]]];
Print["\nUnique pairwise xyz distances: ", udists];

(* 黄金比との比較 *)
phi = N[(1 + Sqrt[5])/2];
Print["φ = ", Round[phi, 0.001], "  φ² = ", Round[phi^2, 0.001],
      "  2/φ = ", Round[2/phi, 0.001]];

(* オービット別に色分けして重心位置を表示 *)
orbit4  = {1, 15, 17, 97};
orbit3  = {4, 18, 51};
singles = {2, 10, 19, 100, 115};

Print["\n=== Orbit positions ==="];
Print["Orbit {1,15,17,97}:"];
Do[pt=cC[[c2]]; Print["  C2=",c2,": xyz=",Round[pt[[{1,2,3}]],0.001],
   " r_xy=",Round[Sqrt[pt[[1]]^2+pt[[2]]^2],0.001],
   " z=",Round[pt[[3]],0.001],
   " θ=",Round[ArcTan[pt[[1]],pt[[2]]]*180/Pi,0.1],"°"],
  {c2, orbit4}];
Print["Orbit {4,18,51}:"];
Do[pt=cC[[c2]]; Print["  C2=",c2,": xyz=",Round[pt[[{1,2,3}]],0.001],
   " r_xy=",Round[Sqrt[pt[[1]]^2+pt[[2]]^2],0.001],
   " z=",Round[pt[[3]],0.001],
   " θ=",Round[ArcTan[pt[[1]],pt[[2]]]*180/Pi,0.1],"°"],
  {c2, orbit3}];
Print["Singletons:"];
Do[pt=cC[[c2]]; Print["  C2=",c2,": xyz=",Round[pt[[{1,2,3}]],0.001],
   " r_xy=",Round[Sqrt[pt[[1]]^2+pt[[2]]^2],0.001],
   " z=",Round[pt[[3]],0.001],
   " θ=",Round[ArcTan[pt[[1]],pt[[2]]]*180/Pi,0.1],"°"],
  {c2, singles}];

(* 3D 散布図 *)
pts3D = {
  {Red,    PointSize[0.04], Point[cC[[#, {1,2,3}]] & /@ orbit4]},
  {Blue,   PointSize[0.04], Point[cC[[#, {1,2,3}]] & /@ orbit3]},
  {Black,  PointSize[0.04], Point[cC[[#, {1,2,3}]] & /@ singles]}
};
gr = Graphics3D[pts3D,
  Axes -> True, AxesLabel -> {"x","y","z"},
  Boxed -> True,
  PlotLabel -> "C2 centroids in xyz (Red={1,15,17,97}, Blue={4,18,51}, Black=singletons)",
  ImageSize -> 400];
Print[gr];
