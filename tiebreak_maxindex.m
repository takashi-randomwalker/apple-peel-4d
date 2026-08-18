(* =====================================================================
   tiebreak_maxindex.m
   On a 300-root sample, max-index was valid at every root of both MIXED
   uniform polytopes, while the queue rule the paper uses managed 1% and
   96%.  Check max-index over ALL 2640 roots of each, so the claim rests
   on a full sweep rather than a sample.

   min-index is run alongside as a second full sweep, since on the
   sample it behaved differently on the two polytopes (99.7% and 80.7%)
   and that difference is worth pinning down exactly.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];
logFile = baseDir <> "tiebreak_maxindex_progress.txt";

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

layersAndOptions[r_] := Module[{dist, q, head, p},
  dist = ConstantArray[-1, nCP]; dist[[r]] = 0;
  q = {r}; head = 1;
  While[head <= Length[q],
   p = q[[head]]; head++;
   Do[If[dist[[nb]] < 0, dist[[nb]] = dist[[p]] + 1; AppendTo[q, nb]],
    {nb, adjListP[[p]]}]];
  {dist, Table[If[i === r, {},
     Select[adjListP[[i]], dist[[#]] === dist[[i]] - 1 &]], {i, nCP}]}];

checkRoot[r_, sel_] := Module[
  {dist, opts, rest, parent, tfsR, unf3D, bMin, bMax, cent, radius,
   cand, ne, hits, treeSet, order},
  {dist, opts} = layersAndOptions[r];
  rest = Complement[Range[nCP], {r}];
  parent = ConstantArray[0, nCP];
  Do[parent[[v]] = sel[opts[[v]]], {v, rest}];
  versP = satFaceUp[r, versCenteredP];
  order = SortBy[rest, dist[[#]] &];
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
  If[Length[cand] === 0, Return[0]];
  ne = Association@Table[
     c -> satNormalsEdges[unf3D[[c]], satFaceTri[[c]], satEdgePair[[c]]],
     {c, Union @@ cand}];
  hits = Select[cand, Function[p,
     satOverlapQ[unf3D[[p[[1]]]], ne[p[[1]]][[1]], ne[p[[1]]][[2]],
                 unf3D[[p[[2]]]], ne[p[[2]]][[1]], ne[p[[2]]][[2]]]]];
  Length[hits]];

runFull[sym_, label_] := Module[{full, out},
  setup[sym];
  full = SelectFirst[h4All, #["sym"] === sym &];
  say["\n", StringRepeat["-", 70]];
  say[label, "  ", sym, "   ", nCP, " cells"];
  say["  queue rule (all roots, from earlier run): ",
      full["valid"], "/", full["cells"]];
  out = Association@Table[
    Module[{t, res, nv},
     {t, res} = AbsoluteTiming@Table[checkRoot[r, sel[[2]]], {r, nCP}];
     nv = Count[res, 0];
     say["  ", StringPadRight[sel[[1]], 10], StringPadLeft[ToString[nv], 6],
         "/", nCP, "  (", Round[100. nv/nCP, 0.01], "%)",
         "   overlaps max ", Max[res], "   ", Round[t, 0.1], " s"];
     sel[[1]] -> <|"valid" -> nv, "overlaps" -> res|>],
    {sel, {{"max-index", Max}, {"min-index", Min}}}];
  <|"sym" -> sym, "label" -> label, "cells" -> nCP, "res" -> out|>];

say[StringRepeat["=", 70]];
say["max-index and min-index over ALL roots of the two MIXED polytopes"];
say[StringRepeat["=", 70]];

tbMax = {
  runFull["x5o3x3x", "runcitruncated 600-cell"],
  runFull["x5x3o3x", "runcitruncated 120-cell"]};

DumpSave[baseDir <> "tiebreak_maxindex.mx", tbMax];
say["\nsaved to tiebreak_maxindex.mx"];
