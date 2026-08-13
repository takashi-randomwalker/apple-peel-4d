(* =====================================================================
   validate_sat.m
   Re-run every root of the 30 A4/B4/F4 polytopes with the SAT test and
   demand the RegionMeasure results back, per root: same valid count and
   the same number of overlapping pairs for each individual root.
   Anything less and the port is wrong.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];

Get[baseDir <> "wythoff_gen_all.mx"];        (* w4All *)
Get[baseDir <> "wythoff_bfs_allroots.mx"];   (* w4BFSAll, RegionMeasure *)
Module[{src, cut},
  src = Import[baseDir <> "face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];
Get[baseDir <> "bfs_sat.m"];

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

Print[StringRepeat["=", 96]];
Print["SAT validation against wythoff_bfs_allroots.mx (RegionMeasure), all roots"];
Print[StringRepeat["=", 96]];
Print[StringPadRight["symbol", 10], StringPadLeft["cells", 6],
      StringPadLeft["valid", 7], StringPadLeft["was", 6],
      "  overlaps/root   ", StringPadLeft["SAT s", 8],
      StringPadLeft["was s", 8], "  verdict"];

res = Table[
  Module[{rec, old, t, out, nValid, sameOv, verdict},
   old = w4BFSAll[[k]];
   rec = SelectFirst[w4All, #["sym"] === old["sym"] &];
   setup[rec];
   {t, out} = AbsoluteTiming[Table[checkRootSAT[r], {r, nCP}]];
   nValid = Count[out, {_, _, _, True}];
   (* per-root overlap counts must agree, not just the totals *)
   sameOv = out[[All, 3]] === old["results"][[All, 3]];
   verdict = If[nValid === old["valid"] && sameOv, "MATCH",
              If[nValid === old["valid"], "counts ok, overlaps DIFFER",
                 "DIFFER"]];
   Print[StringPadRight[old["sym"], 10], StringPadLeft[ToString[nCP], 6],
     StringPadLeft[ToString[nValid], 7], StringPadLeft[ToString[old["valid"]], 6],
     "                   ", StringPadLeft[ToString[Round[t, 0.1]], 8],
     StringPadLeft[ToString[Round[old["time"], 0.1]], 8], "  ", verdict];
   <|"sym" -> old["sym"], "cells" -> nCP, "valid" -> nValid,
     "was" -> old["valid"], "sameOverlaps" -> sameOv,
     "tSAT" -> t, "tOld" -> old["time"],
     "ok" -> (nValid === old["valid"] && sameOv)|>],
  {k, Length[w4BFSAll]}];

Print["\n", StringRepeat["=", 96]];
Print["  matching : ", Count[res, KeyValuePattern["ok" -> True]], " / ",
      Length[res]];
Print["  speedup  : total ", Round[Total[#["tOld"] & /@ res], 0.1], " s -> ",
      Round[Total[#["tSAT"] & /@ res], 0.1], " s   (x",
      Round[Total[#["tOld"] & /@ res]/Total[#["tSAT"] & /@ res], 0.1], ")"];
Do[Print["  MISMATCH: ", r["sym"], "  valid ", r["valid"], " vs ", r["was"],
         "   per-root overlaps equal: ", r["sameOverlaps"]],
   {r, Select[res, ! #["ok"] &]}];
DumpSave[baseDir <> "sat_validation.mx", res];
