(* validate_no_exemption.wls
   Check that the floating-point SAT pipeline (bfs_sat.m), with the
   P-adjacency exemption removed, reproduces the exact results of
   src/certify for the regular polytopes, tree by tree.

   Usage: wolframscript -file validate_no_exemption.wls POLY TREES.json.gz OUT.json [K]
     POLY   f120.m or f600.m (in _4DData/)
     TREES  parent arrays as written by src/certify/export_random_trees.wls
            (per root a list, 0-based, -1 at the root), or "queue" for T_r
     K      use only the first K trees of each root (default: all)
   Writes the overlap count of every tree, in file order. *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
{polyFile, treesFile, outFile} = $ScriptCommandLine[[2 ;; 4]];
kMax = If[Length[$ScriptCommandLine] >= 5, ToExpression[$ScriptCommandLine[[5]]], Infinity];

Module[{src, cut},
  src = Import[baseDir <> "src/bfsnet/face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];
Get[baseDir <> "src/bfsnet/bfs_sat.m"];
(* loadPolytope, layersAndOptions, unfoldWithParents, validQ *)
Module[{src, a, b},
  src = Import[baseDir <> "src/tiebreak/tiebreak_regular.m", "Text"];
  a = First@First@StringPosition[src, "loadPolytope[file_] :="];
  b = First@First@StringPosition[src, "runPolytope[name_, file_] :="];
  ToExpression[StringTake[src, {a, b - 1}]]];

loadPolytope[polyFile];
counts = If[treesFile === "queue",
  Table[checkRootSAT[r][[3]], {r, nCP}],
  Module[{trees = Import[treesFile, {"GZIP", "RawJSON"}]},
   Flatten@Table[
     Module[{dist = layersAndOptions[r][[1]]},
      Table[First@validQ[r, (# + 1) & /@ (trees[[r]][[k]] /. -1 -> 0), dist],
       {k, Min[kMax, Length[trees[[r]]]]}]],
     {r, Length[trees]}]]];
Export[outFile, counts, "RawJSON"];
Print[polyFile, "  ", treesFile, "  trees ", Length[counts],
      "  invalid ", Count[counts, _?Positive], "  overlaps ", Total[counts]];
