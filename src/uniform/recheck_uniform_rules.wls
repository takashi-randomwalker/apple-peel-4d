(* recheck_uniform_rules.wls
   The rule comparison and the random-tree probe for the convex uniform
   4-polytopes, redone with bfs_sat.m after its two fixes of 2026-11-16
   (the earlier runs, src/tiebreak/tiebreak_{all61,maxindex}.m, used the
   buggy test and restricted themselves to polytopes chosen by the buggy
   classification).  All 64 polytopes, every root.

     max      every cell takes the candidate parent of largest index
     min      ... of smallest index
     random   root 1 only, 10 BFS trees drawn uniformly (each cell's parent
              uniform among its candidates), seed 20261116

   Writes data/uniform_recheck/<rule>/<file>.json, per tree:
   {root, overlapping pairs}; skips files already done.

   Usage: wolframscript -file recheck_uniform_rules.wls RULE FILE.m ... *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
Module[{src, cut},
  src = Import[baseDir <> "src/bfsnet/face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];
Get[baseDir <> "src/bfsnet/bfs_sat.m"];
(* setup[] from recheck_uniform_allroots.wls *)
Module[{src, a, b},
  src = Import[baseDir <> "src/uniform/recheck_uniform_allroots.wls", "Text"];
  a = First@First@StringPosition[src, "setup[file_] :="];
  b = First@First@StringPosition[src, "Do[\n Module[{out"];
  ToExpression[StringTake[src, {a, b - 1}]]];
(* layersAndOptions, unfoldWithParents, validQ from tiebreak_regular.m *)
Module[{src, a, b},
  src = Import[baseDir <> "src/tiebreak/tiebreak_regular.m", "Text"];
  a = First@First@StringPosition[src, "(* distance layers from r"];
  b = First@First@StringPosition[src, "runPolytope[name_, file_] :="];
  ToExpression[StringTake[src, {a, b - 1}]]];

rule = $ScriptCommandLine[[2]];
dir = baseDir <> "data/uniform_recheck/" <> rule <> "/";
If[! DirectoryQ[dir], CreateDirectory[dir]];
SeedRandom[20261116];

check[r_, choose_] := Module[{dist, opts, parent},
  {dist, opts} = layersAndOptions[r];
  parent = Table[If[v === r, 0, choose[opts[[v]]]], {v, nCP}];
  {r, First@validQ[r, parent, dist]}];

Do[
 Module[{out = dir <> file <> ".json", t, res},
  If[FileExistsQ[out], Continue[]];
  setup[file];
  {t, res} = AbsoluteTiming@Switch[rule,
     "max", Table[check[r, Max], {r, nCP}],
     "min", Table[check[r, Min], {r, nCP}],
     "random", Table[check[1, RandomChoice], {10}]];
  Export[out, res, "RawJSON", "Compact" -> True];
  Print[DateString[{"Hour", ":", "Minute"}], "  ", rule, "  ", StringPadRight[file, 34],
        StringPadLeft[ToString[nCP], 5], " cells  valid ",
        Count[res, {_, 0}], "/", Length[res], "  (", Round[t], " s)"]],
 {file, $ScriptCommandLine[[3 ;;]]}];
