(* export_sector120_trees.wls
   Write the 2048 D_5-invariant BFS trees of the 120-cell
   (src/symmetry/sector120_nets.wls) as parent arrays, so that
   certify_bfs_exact.py --trees= can check them exactly.

   Output: data/certify/sector_trees_120.json, a one-element list (root 0,
   the paper's cell 1) holding the parent arrays, 0-based with -1 at the
   root.  Run as: certify_bfs_exact.py 120 0 --trees=...

   Usage: wolframscript -file export_sector120_trees.wls *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
Get[baseDir <> "data/sector120_setup.mx"];
img[pm_, x_] := PermutationReplace[x, pm];
nCP = Length[cand];
orbits = DeleteDuplicates[Table[Union[Table[img[h, v], {h, H}]], {v, nCP}]];
reps = First /@ orbits;
fixedCand[v_] := With[{st = Select[H, img[#, v] === v &]},
   Select[cand[[v]], Function[c, AllTrue[st, img[#, c] === c &]]]];
nonRoot = DeleteCases[reps, r]; choices = fixedCand /@ nonRoot;
buildPar[sel_] := Module[{pp = <||>},
  MapThread[Function[{v, c}, Do[pp[img[h, v]] = img[h, c], {h, H}]], {nonRoot, sel}]; pp];
sels = Tuples[choices];
out = Table[Module[{pp = buildPar[s]},
    Table[If[v === r, -1, pp[v] - 1], {v, nCP}]], {s, sels}];
Print["root ", r, "  |H| = ", Length[H], "  trees ", Length[out],
      "  complete ", AllTrue[out, FreeQ[#, _Missing] &]];
Export[baseDir <> "data/certify/sector_trees_120.json", {out}, "RawJSON",
       "Compact" -> True];
