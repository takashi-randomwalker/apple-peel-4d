(* =====================================================================
   export_allroots_csv.m
   One row per root cell, for every convex uniform 4-polytope:
   24,487 rows.  This is the full result set behind Section 4 of the
   paper; the .mx files that hold it are not in the repository.

   Writes uniform_nets/allroots.csv
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];

Get[baseDir <> "wythoff_bfs_allroots.mx"];   (* w4BFSAll : 30, volume test *)
Get[baseDir <> "h4_bfs_allroots.mx"];        (* h4All    : 15, SAT         *)
Get[baseDir <> "special_bfs_allroots.mx"];   (* sp4BFS   :  2, volume test *)
Get[baseDir <> "prisms_bfs_allroots.mx"];    (* pr4BFS   : 17, SAT         *)

tNotation[sym_] := Module[{nodes, marks, rings},
  nodes = StringCases[sym, "x" | "o"];
  marks = StringCases[sym, DigitCharacter];
  If[Length[nodes] =!= 4 || Length[marks] =!= 3, Return[sym]];
  rings = Flatten@Position[nodes, "x"] - 1;
  "t_{" <> StringRiffle[ToString /@ rings, ","] <> "}{" <>
   StringRiffle[marks, ","] <> "}"];

(* the two families that are hard to tell apart, named explicitly *)
knownNames = <|
  "x3o3o3o" -> "5-cell", "o3o3o3x" -> "5-cell",
  "x4o3o3o" -> "tesseract", "o4o3o3x" -> "16-cell",
  "x3o4o3o" -> "24-cell", "o3o4o3x" -> "24-cell",
  "x5o3o3o" -> "120-cell", "o5o3o3x" -> "600-cell",
  "x3o4o3x" -> "runcinated 24-cell",
  "x5o3x3x" -> "runcitruncated 600-cell",
  "x5x3o3x" -> "runcitruncated 120-cell"|>;

rowsFor[recs_, family_, test_] := Join @@ Table[
   Module[{sym = r["sym"], nm},
    nm = Lookup[knownNames, sym, tNotation[sym]];
    Table[
     {sym, nm, family, r["cells"], res[[1]],
      If[res[[4]], "valid", "self-intersecting"], res[[3]], test},
     {res, r["results"]}]],
   {r, recs}];

rows = Join[
  rowsFor[w4BFSAll, "Wythoffian A4/B4/F4", "volume"],
  rowsFor[h4All,    "Wythoffian H4",       "sat"],
  rowsFor[sp4BFS,   "diminished 600-cell", "volume"],
  rowsFor[pr4BFS,   "prism",               "sat"]];

Print["rows: ", Length[rows]];
Print["polytopes: ", Length[DeleteDuplicates[rows[[All, 1]]]]];
Print["valid: ", Count[rows[[All, 6]], "valid"],
      "   self-intersecting: ", Count[rows[[All, 6]], "self-intersecting"]];
Print["by verdict per polytope:"];
Do[Module[{g = GatherBy[Select[rows, #[[1]] === s &], #[[6]] &]},
   If[Length[g] > 1 || g[[1, 1, 6]] =!= "valid",
    Print["   ", StringPadRight[s, 12],
     Count[Select[rows, #[[1]] === s &][[All, 6]], "valid"], "/",
     Length[Select[rows, #[[1]] === s &]], " valid"]]],
 {s, DeleteDuplicates[rows[[All, 1]]]}];

Export[baseDir <> "uniform_nets/allroots.csv",
  Prepend[rows, {"coxeter_symbol", "name", "family", "cells", "root",
                 "verdict", "intersecting_pairs", "intersection_test"}]];

Print["\nwrote uniform_nets/allroots.csv  (",
  Round[FileByteCount[baseDir <> "uniform_nets/allroots.csv"]/1024.^2, 0.01],
  " MB)"];
