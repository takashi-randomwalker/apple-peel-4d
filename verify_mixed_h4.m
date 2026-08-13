(* =====================================================================
   verify_mixed_h4.m
   The two MIXED H4 members are the only results that break the
   all-or-nothing pattern, and x5x3o3x fails with a single overlapping
   pair, which is exactly the regime where the SAT epsilon could be
   deciding the answer.  Re-test a sample of roots with the independent
   RegionMeasure/RegionIntersection criterion.

   Also re-tests the SAT verdict at a few eps values, to see whether the
   failing roots are robust or sit on the threshold.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];

Get[baseDir <> "wythoff_gen_all.mx"];
Get[baseDir <> "h4_bfs_allroots.mx"];
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

mixed = Select[h4All, 0 < #["valid"] < #["cells"] &];

Do[
 Module[{rec, res, bad, good, sample, epsList},
  rec = SelectFirst[w4All, #["sym"] === m["sym"] &];
  setup[rec];
  res  = m["results"];
  bad  = Select[res, ! #[[4]] &][[All, 1]];
  good = Select[res, #[[4]] &][[All, 1]];
  Print["\n", StringRepeat["=", 84]];
  Print[m["sym"], "   ", m["cells"], " cells   SAT: valid ", m["valid"],
        "/", m["cells"], "   overlap tally ", m["overlapTally"]];
  Print[StringRepeat["=", 84]];

  (* 5 failing and 3 passing roots, spread over the range *)
  sample = Join[
    Take[bad, {1, Min[Length[bad], 5 Ceiling[Length[bad]/5]],
               Max[1, Ceiling[Length[bad]/5]]}][[1 ;; UpTo[5]]],
    Take[good, UpTo[3]]];

  Print[StringPadRight["root", 8], StringPadRight["SAT", 22],
        StringPadRight["RegionMeasure", 22], "agree?"];
  Do[
   Module[{s, v},
    s = checkRootSAT[r];
    v = checkRootV2[r];
    Print[StringPadRight[ToString[r], 8],
      StringPadRight[ToString[s[[3]]] <> " overlaps -> " <>
        If[s[[4]], "VALID", "INVALID"], 22],
      StringPadRight[ToString[v[[3]]] <> " overlaps -> " <>
        If[v[[4]], "VALID", "INVALID"], 22],
      If[s[[4]] === v[[4]] && s[[3]] === v[[3]], "yes",
       If[s[[4]] === v[[4]], "verdict only", "*** NO ***"]]]],
   {r, sample}];

  (* epsilon sensitivity on the first failing root *)
  epsList = {10^-4, 10^-5, 10^-6, 10^-7, 10^-8, 10^-10};
  Print["\n  eps sensitivity at root ", First[bad], ":"];
  Do[satEps = e;
   Print["    eps=", ScientificForm[N@e, 2], "   overlaps = ",
         checkRootSAT[First[bad]][[3]]],
   {e, epsList}];
  satEps = 10^-6],
 {m, mixed}];
