(* =====================================================================
   tiebreak_regular.m
   Does Theorem 1 depend on the tie-breaking rule?

   The BFS layering -- the distance from the root -- is canonical: it is
   the same for every tie-breaking.  A BFS tree is therefore exactly a
   choice, for each non-root cell v, of a parent among the neighbors of
   v at distance dist(v)-1.  So the BFS trees rooted at r are in
   bijection with

       Prod_v  |N(v) cap layer(dist(v)-1)|

   and we can both count them exactly and sample them uniformly, without
   reference to any rule.  This replaces "shuffle the queue and hope",
   which does not obviously reach every tree.

   For each root we enumerate every BFS tree when there are few, and
   otherwise draw a uniform sample.  Results go to
   tiebreak_regular.mx and tiebreak_progress.txt.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];
logFile = baseDir <> "tiebreak_progress.txt";
SeedRandom[20260817];

enumCap = 40000;  (* enumerate all BFS trees up to this many *)
nSample = 20;     (* otherwise sample this many *)

Module[{src, cut},
  src = Import[baseDir <> "face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];
Get[baseDir <> "bfs_sat.m"];

say[s__] := (Print[s];
  PutAppend[OutputForm[StringJoin[ToString /@ {s}]], logFile]);

loadPolytope[file_] := Module[{raw},
  raw = Get[baseDir <> "_4DData/" <> file];
  versP  = N[raw[[1]], MachinePrecision];
  facesP = raw[[4]];
  cellsP = raw[[6]];
  nCP    = Length[cellsP];
  versCenteredP = (# - Mean[versP]) & /@ versP;
  adjListP = Table[
    Select[Range[nCP],
      Function[j, j =!= i &&
        Length[Intersection[cellsP[[i]], cellsP[[j]]]] >= 1]], {i, nCP}];
  adjSetP = Association@Flatten[
    Table[Sort[{i, j}] -> True, {i, nCP}, {j, adjListP[[i]]}], 1];
  satPrepCells[]];

(* distance layers from r, and the parent options they induce *)
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

(* unfold with a prescribed parent for every cell *)
unfoldWithParents[r_, parent_, dist_] := Module[{order, tfs},
  order = SortBy[Complement[Range[nCP], {r}], dist[[#]] &];
  tfs = <|r -> idAff4|>;
  Do[tfs[v] = unfoldTFP[parent[[v]], v, tfs[parent[[v]]]], {v, order}];
  tfs];

(* validity of one prescribed tree, via the SAT machinery *)
validQ[r_, parent_, dist_] := Module[
  {tfsR, unf3D, bMin, bMax, cent, radius, cand, ne, hits, treeSet},
  versP = satFaceUp[r, versCenteredP];
  tfsR  = unfoldWithParents[r, parent, dist];
  unf3D = Table[
    Module[{tf = tfsR[i]},
     Take[#, 3] & /@ (versP[[satCellIdx[[i]]]] . Transpose[tf[[1]]] +
                      Threaded[tf[[2]]])], {i, nCP}];
  (* pairs joined by THIS tree are exempt, not the cell-adjacent ones *)
  treeSet = Association@Table[Sort[{v, parent[[v]]}] -> True,
     {v, Complement[Range[nCP], {r}]}];
  bMin = (Min /@ Transpose[#]) & /@ unf3D;
  bMax = (Max /@ Transpose[#]) & /@ unf3D;
  cent = Mean /@ unf3D;
  radius = 2 Max[MapThread[Function[{p, c}, Max[Norm[# - c] & /@ p]],
                           {unf3D, cent}]];
  cand = satBBoxPairs[cent, bMin, bMax, radius];
  cand = Select[cand, ! KeyExistsQ[treeSet, #] &];
  If[Length[cand] === 0, Return[{0, True}]];
  ne = Association@Table[
     c -> satNormalsEdges[unf3D[[c]], satFaceTri[[c]], satEdgePair[[c]]],
     {c, Union @@ cand}];
  hits = Select[cand, Function[p,
     satOverlapQ[unf3D[[p[[1]]]], ne[p[[1]]][[1]], ne[p[[1]]][[2]],
                 unf3D[[p[[2]]]], ne[p[[2]]][[1]], ne[p[[2]]][[2]]]]];
  {Length[hits], Length[hits] === 0}];

runPolytope[name_, file_] := Module[
  {t, res, nTreesAll, checked, bad, mode},
  loadPolytope[file];
  say["\n", StringRepeat["-", 72]];
  say[name, "  (", nCP, " cells)"];
  {t, res} = AbsoluteTiming@Table[
    Module[{dist, opts, nTrees, trees, out},
     {dist, opts} = layersAndOptions[r];
     nTrees = Times @@ (Length /@ Delete[opts, r]);
     If[nTrees <= enumCap,
      mode = "all";
      trees = Tuples[Delete[opts, r]],
      mode = "sample";
      trees = Table[RandomChoice /@ Delete[opts, r], {nSample}]];
     out = Table[
       Module[{parent = ConstantArray[0, nCP], rest},
        rest = Complement[Range[nCP], {r}];
        Do[parent[[rest[[k]]]] = tr[[k]], {k, Length[rest]}];
        validQ[r, parent, dist]],
       {tr, trees}];
     {r, nTrees, Length[out], Count[out, {_, True}], mode}],
    {r, nCP}];
  nTreesAll = res[[All, 2]];
  checked   = Total[res[[All, 3]]];
  bad       = Total[res[[All, 3]] - res[[All, 4]]];
  say["  BFS trees per root: min ", Min[nTreesAll], "  max ", Max[nTreesAll]];
  say["  mode: ", If[MemberQ[res[[All, 5]], "sample"],
                     "sampled " <> ToString[nSample] <> " per root",
                     "EXHAUSTIVE (all BFS trees)"]];
  say["  trees checked: ", checked, "   invalid: ", bad,
      If[bad === 0, "   -> all valid", "   *** SOME INVALID ***"]];
  say["  time ", Round[t, 0.1], " s"];
  <|"name" -> name, "cells" -> nCP, "perRoot" -> res,
    "checked" -> checked, "invalid" -> bad|>];

say[StringRepeat["=", 72]];
say["Tie-break independence: BFS trees other than the index-order one"];
say["enumerate all when <= ", enumCap, " per root, else sample ", nSample];
say[StringRepeat["=", 72]];

tbAll = {
  runPolytope["5-cell",   "f5.m"],
  runPolytope["8-cell",   "f8.m"],
  runPolytope["16-cell",  "f16.m"],
  runPolytope["24-cell",  "f24.m"],
  runPolytope["120-cell", "f120.m"],
  runPolytope["600-cell", "f600.m"]};

DumpSave[baseDir <> "tiebreak_regular.mx", tbAll];
say["\n", StringRepeat["=", 72]];
say["TOTAL trees checked: ", Total[#["checked"] & /@ tbAll],
    "   invalid: ", Total[#["invalid"] & /@ tbAll]];
say["saved to tiebreak_regular.mx"];
