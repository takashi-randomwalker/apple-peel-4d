(* =====================================================================
   tiebreak_all61.m
   The three-way classification of Section 4 was produced with one rule.
   The two MIXED polytopes turned out to move all the way to 100% under
   another rule, so the question is whether the 61 that came out ALL
   VALID are equally rule-dependent, or whether they are robust the way
   the regular polytopes are.

   For each of the 61 (all 64 except x3o4o3x, x5o3x3x, x5x3o3x):
     - every root under max-index
     - every root under min-index
     - root 1 under 10 uniformly random BFS trees, as a cheap probe of
       robustness across the whole family rather than two more rules

   Results are appended to tiebreak_all61.mx after each polytope, so a
   partial run is still usable; polytopes are done smallest first.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];
logFile = baseDir <> "tiebreak_all61_progress.txt";
SeedRandom[20260818];

nRandTrees = 10;

Get[baseDir <> "wythoff_gen_all.mx"];
Get[baseDir <> "uniform4D_prisms.m"];   (* also loads uniform4D_special.m *)
Module[{src, cut},
  src = Import[baseDir <> "face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];
Get[baseDir <> "bfs_sat.m"];

say[s__] := (Print[s];
  PutAppend[OutputForm[StringJoin[ToString /@ {s}]], logFile]);

setupData[data_] := Module[{faceToCells, pairs},
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

nOverlaps[r_, parent_, dist_] := Module[
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
  If[Length[cand] === 0, Return[0]];
  ne = Association@Table[
     c -> satNormalsEdges[unf3D[[c]], satFaceTri[[c]], satEdgePair[[c]]],
     {c, Union @@ cand}];
  hits = Select[cand, Function[p,
     satOverlapQ[unf3D[[p[[1]]]], ne[p[[1]]][[1]], ne[p[[1]]][[2]],
                 unf3D[[p[[2]]]], ne[p[[2]]][[1]], ne[p[[2]]][[2]]]]];
  Length[hits]];

byRule[r_, sel_] := Module[{dist, opts, rest, parent},
  {dist, opts} = layersAndOptions[r];
  rest = Complement[Range[nCP], {r}];
  parent = ConstantArray[0, nCP];
  Do[parent[[v]] = sel[opts[[v]]], {v, rest}];
  nOverlaps[r, parent, dist]];

randomTrees[r_, n_] := Module[{dist, opts, rest},
  {dist, opts} = layersAndOptions[r];
  rest = Complement[Range[nCP], {r}];
  Table[Module[{parent = ConstantArray[0, nCP]},
    Do[parent[[v]] = RandomChoice[opts[[v]]], {v, rest}];
    nOverlaps[r, parent, dist]], {n}]];

skip = {"x3o4o3x", "x5o3x3x", "x5x3o3x"};

jobs = SortBy[
  Join[
   Table[<|"sym" -> r["sym"], "data" -> r["data"], "cells" -> r["counts"][[4]]|>,
     {r, DeleteDuplicatesBy[w4All, {#["counts"], #["hist"]} &]}],
   Table[<|"sym" -> r["sym"], "data" -> r["data"], "cells" -> Length[r["data"][[6]]]|>,
     {r, sp4All}],
   Table[<|"sym" -> StringReplace[r["name"], " " -> "_"], "data" -> r["data"],
           "cells" -> Length[r["data"][[6]]]|>, {r, pr4All}]],
  #["cells"] &];
jobs = Select[jobs, ! MemberQ[skip, #["sym"]] &];

say[StringRepeat["=", 92]];
say["Are the 61 ALL VALID polytopes robust to the tie-breaking rule?"];
say[Length[jobs], " polytopes, ", Total[#["cells"] & /@ jobs], " roots per rule"];
say[StringRepeat["=", 92]];
say[StringPadRight["symbol", 26], StringPadLeft["cells", 6],
    StringPadLeft["max-idx", 12], StringPadLeft["min-idx", 12],
    "   rand@root1"];

all61 = {};
Do[
 Module[{j = jobs[[k]], t, mx, mn, rnd},
  setupData[j["data"]];
  {t, mx} = AbsoluteTiming@Table[byRule[r, Max], {r, nCP}];
  mn = Table[byRule[r, Min], {r, nCP}];
  rnd = randomTrees[1, nRandTrees];
  say[StringPadRight[j["sym"], 26], StringPadLeft[ToString[nCP], 6],
    StringPadLeft[ToString[Count[mx, 0]] <> "/" <> ToString[nCP], 12],
    StringPadLeft[ToString[Count[mn, 0]] <> "/" <> ToString[nCP], 12],
    "   ", Count[rnd, 0], "/", nRandTrees,
    If[Count[mx, 0] < nCP || Count[mn, 0] < nCP, "   *** NOT ALL VALID ***", ""]];
  AppendTo[all61, <|"sym" -> j["sym"], "cells" -> nCP,
     "maxIdx" -> mx, "minIdx" -> mn, "rand" -> rnd|>];
  DumpSave[baseDir <> "tiebreak_all61.mx", all61]],
 {k, Length[jobs]}];

say["\n", StringRepeat["=", 92]];
say["polytopes where max-index is not valid at every root: ",
    Count[all61, x_ /; Count[x["maxIdx"], 0] < x["cells"]]];
say["polytopes where min-index is not valid at every root: ",
    Count[all61, x_ /; Count[x["minIdx"], 0] < x["cells"]]];
say["polytopes where some random tree at root 1 fails: ",
    Count[all61, x_ /; Count[x["rand"], 0] < nRandTrees]];
say["saved to tiebreak_all61.mx"];
