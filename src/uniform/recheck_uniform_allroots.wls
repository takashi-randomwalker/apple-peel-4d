(* recheck_uniform_allroots.wls
   Re-run the face-rotation BFS check of the convex uniform 4-polytopes
   (tree T_r, every root) with bfs_sat.m after its two fixes of 2026-11-16:
   pairs adjacent in P are no longer exempt, and edge directions with a
   zero component are no longer dropped.  Every earlier result in
   uniform_nets/allroots.csv came from the buggy test.

   Reads _4DData_uniform/<file>; writes data/uniform_recheck/<file>.json
   (per root: root, candidate pairs, overlapping pairs) and skips files
   already done, so several copies can share the list and a crash loses
   at most one polytope.

   Usage: wolframscript -file recheck_uniform_allroots.wls FILE.m ... *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
Module[{src, cut},
  src = Import[baseDir <> "src/bfsnet/face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];
Get[baseDir <> "src/bfsnet/bfs_sat.m"];

(* same globals as setup[] in run_h4_allroots.m; adjacency rebuilt from
   the faces, ascending, so T_r is the tree of the original runs *)
setup[file_] := Module[{raw, faceToCells, pairs},
  raw    = Get[baseDir <> "_4DData_uniform/" <> file];
  versP  = N[raw[[1]], MachinePrecision];
  facesP = raw[[3]];
  cellsP = raw[[4]];
  nCP    = Length[cellsP];
  versCenteredP = (# - Mean[versP]) & /@ versP;
  faceToCells = GroupBy[
    Flatten@MapIndexed[Function[{cf, k}, {# -> First[k]} & /@ cf], cellsP],
    First -> Last];
  pairs = Select[Values[faceToCells], Length[#] === 2 &];
  adjListP = Union /@ Lookup[
    GroupBy[Join[pairs, Reverse /@ pairs], First -> Last], Range[nCP], {}];
  adjSetP = Association@Flatten[
    Table[Sort[{i, j}] -> True, {i, nCP}, {j, adjListP[[i]]}], 1];
  If[Sort /@ adjListP =!= Sort /@ raw[[5]], Print["  cellAdj of the file disagrees"]];
  satPrepCells[]];

Do[
 Module[{out = baseDir <> "data/uniform_recheck/" <> file <> ".json", t, res},
  If[FileExistsQ[out], Continue[]];
  setup[file];
  {t, res} = AbsoluteTiming[Table[Most[checkRootSAT[r]], {r, nCP}]];
  Export[out, res, "RawJSON", "Compact" -> True];
  Print[DateString[{"Hour", ":", "Minute"}], "  ", StringPadRight[file, 34],
        StringPadLeft[ToString[nCP], 5], " cells  valid ",
        Count[res, {_, _, 0}], "/", nCP, "  (", Round[t], " s)"]],
 {file, $ScriptCommandLine[[2 ;;]]}];
