(* =====================================================================
   check_runcinated24.m
   The runcinated 24-cell x3o4o3x was the only Wythoffian member of the
   A4/B4/F4 families whose face-rotation BFS net self-intersects at
   root 1.  Run every root and record how badly it fails.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];

Get[baseDir <> "uniform4D_wythoff.m"];
Module[{src, cut},
  src = Import[baseDir <> "face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];

r24    = w4Find["x3o4o3x"];
data   = r24["data"];
versP  = N[data[[1]], MachinePrecision];
facesP = data[[4]];
cellsP = data[[6]];
nCP    = Length[cellsP];
versCenteredP = (# - Mean[versP]) & /@ versP;

adjListP = Table[
  Select[Range[nCP],
    Function[j, j =!= i &&
      Length[Intersection[cellsP[[i]], cellsP[[j]]]] >= 1]], {i, nCP}];
adjSetP = Association@Flatten[
  Table[Sort[{i, j}] -> True, {i, nCP}, {j, adjListP[[i]]}], 1];

Print["runcinated 24-cell  x3o4o3x"];
Print["  V,E,F,C = ", r24["counts"], "   cells = ", r24["hist"]];
Print["  cell degree: ", Min[Length /@ adjListP], "-", Max[Length /@ adjListP],
      "   degree tally = ", Sort@Tally[Length /@ adjListP]];
Print["  running ", nCP, " roots ..."];

{tAll, res} = AbsoluteTiming[Table[checkRootV2[r], {r, nCP}]];
valid = Select[res, #[[4]] &];

Print["\n  valid roots : ", Length[valid], " / ", nCP];
Print["  overlapping pairs per root: min=", Min[res[[All, 3]]],
      "  max=", Max[res[[All, 3]]],
      "  mean=", Round[Mean[N@res[[All, 3]]], 0.01]];
Print["  bbox candidate pairs per root: min=", Min[res[[All, 2]]],
      "  max=", Max[res[[All, 2]]]];
Print["  tally of true-overlap counts: ", Sort@Tally[res[[All, 3]]]];
If[Length[valid] > 0, Print["  valid roots: ", valid[[All, 1]]]];
Print["  total time: ", Round[tAll, 0.1], "s"];
