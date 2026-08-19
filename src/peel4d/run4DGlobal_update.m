(* ================================================================
   run4DGlobal_update.m
   ================================================================
   現行 peeling4Df4.m（k=2 xy クロス積タイブレーク版）で
   5胞体・8胞体・16胞体・24胞体を再計算する．
   120胞体は ans4DGlobal.mx の既存データを流用する（600胞体は除外）．
   results4DGlobal として ans4DGlobal_v4.mx に保存する．

   使い方 (Mathematica Notebook):
     SetDirectory["/path/to/260324Peeling4D"];
     Get["run4DGlobal_update.m"]
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];
Get[scriptDir <> "peeling4Df4.m"];
dataDir = scriptDir <> "_4DData/";

polytopeNames = <|5->"5-cell", 8->"8-cell", 16->"16-cell",
                  24->"24-cell", 120->"120-cell"|>;

loadPoly[n_] := Module[{raw},
  raw = Get[dataDir <> "f" <> ToString[n] <> ".m"];
  <|"vers"->raw[[1]], "faces"->raw[[4]], "cells"->raw[[6]]|>
];

(* 既存データ読み込み（120胞体を流用）*)
Print["Loading ans4DGlobal.mx (120-cell data preserved) ..."];
Get[scriptDir <> "ans4DGlobal.mx"];  (* defines results4DGlobal *)

Print[""];
Print[StringRepeat["=", 65]];
Print["  Apple-Peel 4D: 5/8/16/24-cell  (k=2 xy tie-break)"];
Print[StringRepeat["=", 65]];

Do[
  poly   = loadPoly[n];
  vers   = poly["vers"]; faces = poly["faces"]; cells = poly["cells"];
  nc     = Length[cells];
  Print["\nProcessing ", polytopeNames[n], " (", nc, " cells) ..."];
  t0     = AbsoluteTime[];

  ans = Table[peeling4Df4[N[vers], faces, cells, top], {top, nc}];

  dt         = Round[AbsoluteTime[] - t0, 0.1];
  allBools   = Flatten[Map[Last, ans, {2}]];
  trueCount  = Count[allBools, True];
  falseCount = Count[allBools, False];
  uniqOrders = DeleteDuplicates[
    Select[Flatten[ans, 1], Last[#]&][[All, 2]]];
  classif    = Which[falseCount == 0, "Perfect",
                     trueCount  == 0, "Impossible",
                     True,            "Possible"];

  Print["  ", dt, " s  |  ", classif,
        "  (True: ", trueCount,
        ", False: ", falseCount,
        ", Unique: ", Length[uniqOrders], ")"];

  (* 旧版との比較 *)
  If[KeyExistsQ[results4DGlobal, n],
    old = results4DGlobal[n];
    Print["  old: ", old["classification"],
          "  (True: ", old["trueCount"],
          ", Unique: ", old["uniqueCount"], ")"]
  ];

  results4DGlobal[n] = <|
    "name"           -> polytopeNames[n],
    "n"              -> n,
    "classification" -> classif,
    "trueCount"      -> trueCount,
    "falseCount"     -> falseCount,
    "uniqueCount"    -> Length[uniqOrders],
    "results"        -> ans
  |>,
  {n, {5, 8, 16, 24}}
];

(* 保存 *)
outFile = scriptDir <> "ans4DGlobal_v4.mx";
DumpSave[outFile, results4DGlobal];
Print["\nSaved: ", outFile];

(* サマリ *)
Print[""];
Print[StringRepeat["=", 65]];
Print["  Summary"];
Print[StringRepeat["=", 65]];
Print[StringPadLeft["Polytope", 10], "  ",
      StringPadLeft["Class",    10], "  ",
      StringPadLeft["True",      6], "  ",
      StringPadLeft["Unique",    6]];
Print[StringRepeat["-", 40]];
Do[
  r = results4DGlobal[n];
  Print[StringPadLeft[r["name"],           10], "  ",
        StringPadLeft[r["classification"], 10], "  ",
        StringPadLeft[ToString[r["trueCount"]],  6], "  ",
        StringPadLeft[ToString[r["uniqueCount"]], 6]],
  {n, {5, 8, 16, 24, 120}}
];

Print["\nDone."];
