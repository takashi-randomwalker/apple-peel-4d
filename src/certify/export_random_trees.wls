(* export_random_trees.wls
   Redraw the random BFS trees of the paper (20 per root for the 120- and
   600-cell, src/tiebreak/tiebreak_regular.m) so that certify_bfs_exact.py
   can check the very same trees exactly.

   Only the draws are reproduced: same seed, same adjacency lists, same
   layersAndOptions, same order of RandomChoice calls.  The 5- to 24-cell
   were enumerated there and consumed no random numbers.

   Writes data/certify/random_trees_{120,600}.json: for each root a list
   of 20 parent arrays, 0-based, with -1 at the root.

   Usage: wolframscript -file export_random_trees.wls *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SeedRandom[20260817];
nSample = 20;

loadPolytope[file_] := Module[{raw},
  raw = Get[baseDir <> "_4DData/" <> file];
  cellsP = raw[[6]];
  nCP    = Length[cellsP];
  adjListP = Table[
    Select[Range[nCP],
      Function[j, j =!= i &&
        Length[Intersection[cellsP[[i]], cellsP[[j]]]] >= 1]], {i, nCP}]];

layersAndOptions[r_] := Module[{dist, q, head, p, opts},
  dist = ConstantArray[-1, nCP]; dist[[r]] = 0;
  q = {r}; head = 1;
  While[head <= Length[q],
   p = q[[head]]; head++;
   Do[If[dist[[nb]] < 0, dist[[nb]] = dist[[p]] + 1; AppendTo[q, nb]],
    {nb, adjListP[[p]]}]];
  opts = Table[
    If[i === r, {},
     Select[adjListP[[i]], dist[[#]] === dist[[i]] - 1 &]], {i, nCP}];
  {dist, opts}];

Do[
  loadPolytope["f" <> name <> ".m"];
  out = Table[
    Module[{opts, trees},
     opts = Last[layersAndOptions[r]];
     trees = Table[RandomChoice /@ Delete[opts, r], {nSample}];
     Table[
      Module[{parent = ConstantArray[0, nCP], rest},
       rest = Complement[Range[nCP], {r}];
       Do[parent[[rest[[k]]]] = tr[[k]], {k, Length[rest]}];
       parent - 1],
      {tr, trees}]],
    {r, nCP}];
  Export[baseDir <> "data/certify/random_trees_" <> name <> ".json", out, "RawJSON",
         "Compact" -> True];
  Print[name, "-cell: ", nCP, " roots x ", nSample, " trees written"],
  {name, {"120", "600"}}];
