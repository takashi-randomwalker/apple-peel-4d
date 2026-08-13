(* =====================================================================
   run_h4_allroots.m
   Face-rotation BFS net check over EVERY root for all 15 Wythoffian
   members of the H4 family, using the SAT test (bfs_sat.m).

   14,040 roots in total.  Estimated ~6-7 h; results are appended to
   h4_bfs_allroots.mx after each polytope so a crash loses at most one.
   Progress goes to h4_allroots_progress.txt as well as stdout.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];
logFile = baseDir <> "h4_allroots_progress.txt";

Get[baseDir <> "wythoff_gen_all.mx"];
Module[{src, cut},
  src = Import[baseDir <> "face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];
Get[baseDir <> "bfs_sat.m"];

say[s__] := (Print[s];
  PutAppend[OutputForm[StringJoin[ToString /@ {s}]], logFile]);

setup[rec_] := Module[{data, faceToCells, pairs},
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
  satPrepCells[]];

h4 = SortBy[Select[w4All, StringContainsQ[#["sym"], "5"] &],
            #["counts"][[4]] &];

say[StringRepeat["=", 96]];
say["H4 family, face-rotation BFS nets, ALL roots (SAT).  ",
    Length[h4], " polytopes, ", Total[#["counts"][[4]] & /@ h4], " roots"];
say[StringRepeat["=", 96]];

h4All = {};
Do[
 Module[{rec = h4[[k]], t, res, nValid, degs},
  setup[rec];
  degs = Sort@Tally[Length /@ adjListP];
  {t, res} = AbsoluteTiming[Table[checkRootSAT[r], {r, nCP}]];
  nValid = Count[res, {_, _, _, True}];
  say[StringPadRight[rec["sym"], 10], StringPadLeft[ToString[nCP], 6],
   " cells   valid ", StringPadLeft[ToString[nValid], 6], "/",
   StringPadRight[ToString[nCP], 6],
   StringPadRight[
    Which[nValid === nCP, "ALL VALID", nValid === 0, "ALL INVALID",
          True, "MIXED"], 13],
   "overlaps min/max ", Min[res[[All, 3]]], "/", Max[res[[All, 3]]],
   "   deg ", degs, "   ", Round[t, 0.1], "s"];
  AppendTo[h4All,
   <|"sym" -> rec["sym"], "counts" -> rec["counts"], "hist" -> rec["hist"],
     "cells" -> nCP, "valid" -> nValid, "degs" -> degs,
     "overlapTally" -> Sort@Tally[res[[All, 3]]],
     "results" -> res, "time" -> t|>];
  DumpSave[baseDir <> "h4_bfs_allroots.mx", h4All]],
 {k, Length[h4]}];

say[""];
say[StringRepeat["=", 96]];
say["SUMMARY"];
say["  polytopes    : ", Length[h4All]];
say["  roots        : ", Total[#["cells"] & /@ h4All]];
say["  valid roots  : ", Total[#["valid"] & /@ h4All]];
say["  ALL VALID    : ", Count[h4All, r_ /; r["valid"] === r["cells"]]];
say["  MIXED        : ", Count[h4All, r_ /; 0 < r["valid"] < r["cells"]]];
say["  ALL INVALID  : ", Count[h4All, r_ /; r["valid"] === 0]];
say["  not ALL VALID: ",
    Select[h4All, #["valid"] =!= #["cells"] &][[All, "sym"]]];
say["  saved to h4_bfs_allroots.mx"];
