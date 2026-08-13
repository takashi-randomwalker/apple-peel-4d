(* =====================================================================
   run_special_bfs.m
   Face-rotation BFS net check over EVERY root for the snub 24-cell and
   the grand antiprism (uniform4D_special.m).  C = 144 and 320, so all
   roots is affordable here.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];

Get[baseDir <> "uniform4D_special.m"];
Module[{src, cut},
  src = Import[baseDir <> "face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];

runAll[rec_] := Module[{data, res, t, faceToCells, pairs, nValid},
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

  {t, res} = AbsoluteTiming[Table[checkRootV2[r], {r, nCP}]];
  nValid = Count[res, {_, _, _, True}];

  Print[StringPadRight[rec["name"], 18],
   StringPadLeft[ToString[nCP], 5], " cells   valid ",
   StringPadLeft[ToString[nValid], 4], "/", StringPadRight[ToString[nCP], 5],
   StringPadRight[
    Which[nValid === nCP, "ALL VALID", nValid === 0, "ALL INVALID",
          True, "MIXED"], 13],
   "overlaps min/max ", Min[res[[All, 3]]], "/", Max[res[[All, 3]]],
   "   cell deg ", Sort@Tally[Length /@ adjListP],
   "   ", Round[t, 0.1], "s"];
  Join[rec, <|"valid" -> nValid, "cells" -> nCP,
              "overlapTally" -> Sort@Tally[res[[All, 3]]],
              "results" -> res|>]];

Print["\n", StringRepeat["=", 100]];
Print["Face-rotation BFS nets, ALL roots: snub 24-cell and grand antiprism"];
Print[StringRepeat["=", 100]];

sp4BFS = runAll /@ sp4All;
DumpSave[baseDir <> "special_bfs_allroots.mx", sp4BFS];
Print["\n  saved to special_bfs_allroots.mx"];
