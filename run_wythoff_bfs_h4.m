(* =====================================================================
   run_wythoff_bfs_h4.m
   1. Cross-check the generated 120-cell and 600-cell against
      _4DData/f120.m and f600.m.
   2. Face-rotation BFS net check, ROOT 1 ONLY, for the H4 family.

   All roots is out of reach here: the 240-cell x3o4o3x needed 229 s for
   240 roots, and the per-root cost grows like C^2, so a 2640-cell member
   would need days.  Root 1 is reported as what it is - a single root.
   In the 30 A4/B4/F4 polytopes the outcome was always all-or-nothing,
   which makes root 1 a strong indicator, but that dichotomy is an
   observation, not a theorem.

   Cell adjacency is built from the ridge incidence (every 2-face lies in
   exactly two cells) instead of the O(C^2) pairwise Intersection scan
   used by face_rotation_net_all4D_v2.m; the two agree by construction.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];

Get[baseDir <> "wythoff_gen_all.mx"];
Module[{src, cut},
  src = Import[baseDir <> "face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];

w4Find[sym_] := SelectFirst[w4All,
   #["sym"] === sym || #["sym"] === StringReverse[sym] &];

(* --- 1. cross-check the two regular H4 members ------------------- *)
Print[StringRepeat["=", 78]];
Print["Cross-check against _4DData/ (V, F, C from parts 1, 4, 6)"];
Print[StringRepeat["=", 78]];
Do[
 Module[{raw = Get[baseDir <> "_4DData/" <> p[[1]]], want, got, r},
  want = Length /@ raw[[{1, 4, 6}]];
  r    = w4Find[p[[2]]];
  got  = r["counts"][[{1, 3, 4}]];
  Print["  ", StringPadRight[p[[1]], 8], StringPadRight[p[[2]], 9],
        " file=", want, " wythoff=", got,
        If[want === got, "   MATCH", "   DIFFER"],
        "   edges ", Length[Union[Sort /@ raw[[2]]]], "/",
        r["counts"][[2]]]],
 {p, {{"f120.m", "x5o3o3o"}, {"f600.m", "o5o3o3x"}}}];

(* --- 2. BFS, root 1 ---------------------------------------------- *)
runRoot1[rec_, cap_] := Module[
  {data, res, t, faceToCells, pairs},
  data   = rec["data"];
  versP  = N[data[[1]], MachinePrecision];
  facesP = data[[4]];
  cellsP = data[[6]];
  nCP    = Length[cellsP];
  versCenteredP = (# - Mean[versP]) & /@ versP;

  (* adjacency from ridge incidence: each 2-face joins exactly two cells *)
  faceToCells = GroupBy[
    Flatten@MapIndexed[Function[{cf, k}, {# -> First[k]} & /@ cf], cellsP],
    First -> Last];
  pairs = Select[Values[faceToCells], Length[#] === 2 &];
  adjListP = Lookup[
    GroupBy[Join[pairs, Reverse /@ pairs], First -> Last],
    Range[nCP], {}];
  adjListP = Union /@ adjListP;
  adjSetP = Association@Flatten[
    Table[Sort[{i, j}] -> True, {i, nCP}, {j, adjListP[[i]]}], 1];

  {t, res} = AbsoluteTiming@TimeConstrained[checkRootV2[1], cap, $Aborted];

  If[res === $Aborted,
   Print[StringPadRight[rec["sym"], 10], StringPadLeft[ToString[nCP], 6],
     " cells   TIMEOUT after ", Round[t], "s"];
   Return[<|"sym" -> rec["sym"], "cells" -> nCP, "status" -> "timeout"|>]];

  Print[StringPadRight[rec["sym"], 10], StringPadLeft[ToString[nCP], 6],
   " cells   ", If[res[[4]], "VALID  ", "INVALID"],
   "   bbox pairs ", StringPadLeft[ToString[res[[2]]], 7],
   "   overlaps ", StringPadLeft[ToString[res[[3]]], 5],
   "   deg ", Min[Length /@ adjListP], "-", Max[Length /@ adjListP],
   "   ", Round[t, 0.1], "s"];
  <|"sym" -> rec["sym"], "cells" -> nCP, "valid" -> res[[4]],
    "bbox" -> res[[2]], "overlaps" -> res[[3]], "status" -> "done"|>];

h4 = SortBy[Select[w4All, StringContainsQ[#["sym"], "5"] &],
            #["counts"][[4]] &];

Print["\n", StringRepeat["=", 78]];
Print["H4 family, face-rotation BFS net at ROOT 1 (", Length[h4], " members)"];
Print[StringRepeat["=", 78]];

h4res = runRoot1[#, 3000] & /@ h4;

Print["\n  done ", Count[h4res, KeyValuePattern["status" -> "done"]], "/",
      Length[h4res],
      "   valid ", Count[h4res, KeyValuePattern["valid" -> True]],
      "   invalid ", Count[h4res, KeyValuePattern["valid" -> False]],
      "   timeout ", Count[h4res, KeyValuePattern["status" -> "timeout"]]];
Print["  invalid at root 1: ",
  Select[h4res, #["status"] === "done" && ! #["valid"] &][[All, "sym"]]];

DumpSave[baseDir <> "wythoff_bfs_h4_root1.mx", h4res];
Print["\n  saved to wythoff_bfs_h4_root1.mx"];
