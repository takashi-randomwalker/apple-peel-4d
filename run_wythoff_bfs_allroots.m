(* =====================================================================
   run_wythoff_bfs_allroots.m
   Face-rotation BFS net check over EVERY root cell, for all Wythoffian
   members of the A4, B4 and F4 families (30 distinct polytopes,
   1893 roots in total).

   Companion to run_wythoff_bfs.m, which only did root 1 (plus all roots
   for the four regular members).  Results are saved to
   wythoff_bfs_allroots.mx as a list of associations.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];

Get[baseDir <> "uniform4D_wythoff.m"];
Module[{src, cut},
  src = Import[baseDir <> "face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];

w4Distinct = SortBy[
  DeleteDuplicatesBy[w4All, {#["counts"], #["hist"]} &],
  #["counts"][[4]] &];

runAllRoots[rec_] := Module[
  {data, res, t, nValid, degs},
  data   = rec["data"];
  versP  = N[data[[1]], MachinePrecision];
  facesP = data[[4]];
  cellsP = data[[6]];
  nCP    = Length[cellsP];
  versCenteredP = (# - Mean[versP]) & /@ versP;

  adjListP = Table[
    Select[Range[nCP],
      Function[j, j =!= i &&
        Length[Intersection[cellsP[[i]], cellsP[[j]]]] >= 1]], {i, nCP}];
  adjSetP = Association@Flatten[
    Table[Sort[{i, j}] -> True, {i, nCP}, {j, adjListP[[i]]}], 1];
  degs = Sort@Tally[Length /@ adjListP];

  {t, res} = AbsoluteTiming[Table[checkRootV2[r], {r, nCP}]];
  nValid = Count[res, {_, _, _, True}];

  Print[
   StringPadRight[rec["sym"], 10],
   StringPadLeft[ToString[nCP], 6], " cells   valid ",
   StringPadLeft[ToString[nValid], 5], "/", StringPadRight[ToString[nCP], 5],
   StringPadRight[
     Which[nValid === nCP, "ALL VALID", nValid === 0, "ALL INVALID",
           True, "MIXED"], 13],
   "overlaps min/max ", StringPadLeft[ToString[Min[res[[All, 3]]]], 4], "/",
   StringPadRight[ToString[Max[res[[All, 3]]]], 5],
   "deg ", StringPadRight[ToString[degs], 22],
   StringPadLeft[ToString[Round[t, 0.1]], 8], "s"];

  <|"sym" -> rec["sym"], "counts" -> rec["counts"], "hist" -> rec["hist"],
    "cells" -> nCP, "valid" -> nValid, "degs" -> degs,
    "overlapMin" -> Min[res[[All, 3]]], "overlapMax" -> Max[res[[All, 3]]],
    "overlapTally" -> Sort@Tally[res[[All, 3]]],
    "results" -> res, "time" -> t|>];

Print[StringRepeat["=", 108]];
Print["Face-rotation BFS nets, ALL roots, A4 + B4 + F4 Wythoffian members"];
Print["  ", Length[w4Distinct], " polytopes, ",
      Total[#["counts"][[4]] & /@ w4Distinct], " roots"];
Print[StringRepeat["=", 108]];

w4BFSAll = runAllRoots /@ w4Distinct;

DumpSave[baseDir <> "wythoff_bfs_allroots.mx", w4BFSAll];

Print["\n", StringRepeat["=", 108]];
Print["SUMMARY"];
Print[StringRepeat["=", 108]];
Print["  polytopes          : ", Length[w4BFSAll]];
Print["  roots              : ", Total[#["cells"] & /@ w4BFSAll]];
Print["  valid roots        : ", Total[#["valid"] & /@ w4BFSAll]];
Print["  ALL VALID polytopes: ",
      Count[w4BFSAll, r_ /; r["valid"] === r["cells"]]];
Print["  MIXED polytopes    : ",
      Count[w4BFSAll, r_ /; 0 < r["valid"] < r["cells"]]];
Print["  ALL INVALID        : ", Count[w4BFSAll, r_ /; r["valid"] === 0]];
Print["\n  not ALL VALID:"];
Do[Print["    ", StringPadRight[r["sym"], 10], r["valid"], "/", r["cells"],
         "   cells=", r["hist"], "   deg=", r["degs"]],
   {r, Select[w4BFSAll, #["valid"] =!= #["cells"] &]}];
Print["\n  saved to wythoff_bfs_allroots.mx"];
