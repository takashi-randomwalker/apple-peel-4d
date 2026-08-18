(* =====================================================================
   tiebreak_mixed.m
   Is the root-dependence of the two MIXED uniform 4-polytopes an
   artifact of the tie-breaking rule?

   x5o3x3x (runcitruncated 600-cell) is valid from 47 of its 2640 roots
   under the index-order rule; x5x3o3x (runcitruncated 120-cell) from
   2570 of 2640.  Both have 2640 cells.  For a sample of roots of each
   kind -- ones the fixed rule calls valid and ones it calls invalid --
   we draw BFS trees uniformly at random (same construction as
   tiebreak_regular.m) and ask what fraction of them are valid.

   Three outcomes are worth distinguishing:
     - a root the rule calls invalid is invalid for every tree we draw
       => root-dependence is a property of the root, not the rule
     - such a root is valid for some trees
       => the rule, not the root, is deciding, and the counts 47/2640
          and 2570/2640 say less than they appear to
     - a root the rule calls valid fails for some tree
       => same conclusion from the other side
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];
logFile = baseDir <> "tiebreak_mixed_progress.txt";
SeedRandom[20260817];

nRoots  = 12;   (* roots sampled from each class *)
nTrees  = 25;   (* BFS trees drawn per root *)

Get[baseDir <> "wythoff_gen_all.mx"];
Get[baseDir <> "h4_bfs_allroots.mx"];
Module[{src, cut},
  src = Import[baseDir <> "face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];
Get[baseDir <> "bfs_sat.m"];

say[s__] := (Print[s];
  PutAppend[OutputForm[StringJoin[ToString /@ {s}]], logFile]);

setup[sym_] := Module[{rec, data, faceToCells, pairs},
  rec = SelectFirst[w4All, #["sym"] === sym &];
  data = rec["data"];
  versP  = N[data[[1]], MachinePrecision];
  facesP = data[[4]];
  cellsP = data[[6]];
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
  satPrepCells[]];

layersAndOptions[r_] := Module[{dist, q, head, p, opts},
  dist = ConstantArray[-1, nCP]; dist[[r]] = 0;
  q = {r}; head = 1;
  While[head <= Length[q],
   p = q[[head]]; head++;
   Do[If[dist[[nb]] < 0, dist[[nb]] = dist[[p]] + 1; AppendTo[q, nb]],
    {nb, adjListP[[p]]}]];
  {dist, Table[If[i === r, {},
     Select[adjListP[[i]], dist[[#]] === dist[[i]] - 1 &]], {i, nCP}]}];

validQ[r_, parent_, dist_] := Module[
  {tfsR, unf3D, bMin, bMax, cent, radius, cand, ne, hits, treeSet, order},
  versP = satFaceUp[r, versCenteredP];
  order = SortBy[Complement[Range[nCP], {r}], dist[[#]] &];
  tfsR = <|r -> idAff4|>;
  Do[tfsR[v] = unfoldTFP[parent[[v]], v, tfsR[parent[[v]]]], {v, order}];
  unf3D = Table[
    Module[{tf = tfsR[i]},
     Take[#, 3] & /@ (versP[[satCellIdx[[i]]]] . Transpose[tf[[1]]] +
                      Threaded[tf[[2]]])], {i, nCP}];
  treeSet = Association@Table[Sort[{v, parent[[v]]}] -> True, {v, order}];
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

sampleRoot[r_] := Module[{dist, opts, rest, out},
  {dist, opts} = layersAndOptions[r];
  rest = Complement[Range[nCP], {r}];
  out = Table[
    Module[{parent = ConstantArray[0, nCP]},
     Do[parent[[v]] = RandomChoice[opts[[v]]], {v, rest}];
     validQ[r, parent, dist]],
    {nTrees}];
  {Count[out, {_, True}], Length[out], Min[out[[All, 1]]], Max[out[[All, 1]]]}];

runMixed[sym_, label_] := Module[{rec, good, bad, pick, res, t},
  setup[sym];
  rec = SelectFirst[h4All, #["sym"] === sym &];
  good = Select[rec["results"], #[[4]] &][[All, 1]];
  bad  = Select[rec["results"], ! #[[4]] &][[All, 1]];
  say["\n", StringRepeat["-", 74]];
  say[label, "  ", sym, "   ", nCP, " cells   fixed rule: ",
      Length[good], " valid / ", Length[bad], " invalid roots"];
  pick = Join[
    {"rule-valid", #} & /@ RandomSample[good, Min[nRoots, Length[good]]],
    {"rule-invalid", #} & /@ RandomSample[bad, Min[nRoots, Length[bad]]]];
  {t, res} = AbsoluteTiming@Table[
    Module[{s = sampleRoot[p[[2]]]},
     say["  ", StringPadRight[p[[1]], 13], " root ",
         StringPadLeft[ToString[p[[2]]], 5], ":  ",
         s[[1]], "/", s[[2]], " random trees valid",
         "   overlaps ", s[[3]], "-", s[[4]]];
     {p[[1]], p[[2]], s}],
    {p, pick}];
  say["  time ", Round[t, 0.1], " s"];
  <|"sym" -> sym, "label" -> label, "res" -> res|>];

say[StringRepeat["=", 74]];
say["Do the MIXED uniform polytopes stay root-dependent under other"];
say["BFS tie-breakings?   ", nRoots, " roots per class, ",
    nTrees, " random trees per root"];
say[StringRepeat["=", 74]];

tbMixed = {
  runMixed["x5o3x3x", "runcitruncated 600-cell"],
  runMixed["x5x3o3x", "runcitruncated 120-cell"]};

DumpSave[baseDir <> "tiebreak_mixed.mx", tbMixed];
say["\nsaved to tiebreak_mixed.mx"];
