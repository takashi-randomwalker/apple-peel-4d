(* =====================================================================
   run_prisms_bfs.m
   Face-rotation BFS net check over EVERY root for the 17 prismatic
   convex uniform 4-polytopes (uniform4D_prisms.m).  530 roots in all.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];

Get[baseDir <> "uniform4D_prisms.m"];
Module[{src, cut},
  src = Import[baseDir <> "face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];
Get[baseDir <> "bfs_sat.m"];

runAll[rec_] := Module[{data, res, t, faceToCells, pairs, nValid, degs},
  data   = rec["data"];
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
  satPrepCells[];
  degs = Sort@Tally[Length /@ adjListP];

  {t, res} = AbsoluteTiming[Table[checkRootSAT[r], {r, nCP}]];
  nValid = Count[res, {_, _, _, True}];

  Print[StringPadRight[rec["name"], 30],
   StringPadLeft[ToString[nCP], 5], " cells   valid ",
   StringPadLeft[ToString[nValid], 4], "/", StringPadRight[ToString[nCP], 5],
   StringPadRight[
    Which[nValid === nCP, "ALL VALID", nValid === 0, "ALL INVALID",
          True, "MIXED"], 13],
   "overlaps min/max ", Min[res[[All, 3]]], "/", Max[res[[All, 3]]],
   "   deg ", degs, "   ", Round[t, 0.1], "s"];
  Join[rec, <|"valid" -> nValid, "cells" -> nCP,
              "degs" -> degs,
              "overlapTally" -> Sort@Tally[res[[All, 3]]],
              "results" -> res, "time" -> t|>]];

Print["\n", StringRepeat["=", 108]];
Print["Face-rotation BFS nets, ALL roots: the 17 prismatic uniform 4-polytopes (",
      Total[#["counts"][[4]] & /@ pr4All], " roots)"];
Print[StringRepeat["=", 108]];

pr4BFS = runAll /@ pr4All;
DumpSave[baseDir <> "prisms_bfs_allroots.mx", pr4BFS];

Print["\n", StringRepeat["=", 108]];
Print["  polytopes    : ", Length[pr4BFS]];
Print["  roots        : ", Total[#["cells"] & /@ pr4BFS]];
Print["  valid roots  : ", Total[#["valid"] & /@ pr4BFS]];
Print["  ALL VALID    : ", Count[pr4BFS, r_ /; r["valid"] === r["cells"]]];
Print["  MIXED        : ", Count[pr4BFS, r_ /; 0 < r["valid"] < r["cells"]]];
Print["  ALL INVALID  : ", Count[pr4BFS, r_ /; r["valid"] === 0]];
Print["  not ALL VALID: ",
  Column[{#["name"], #["valid"], "/", #["cells"]} & /@
    Select[pr4BFS, #["valid"] =!= #["cells"] &]]];
Print["  saved to prisms_bfs_allroots.mx"];
