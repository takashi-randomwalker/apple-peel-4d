(* ================================================================
   run4DGlobal_tiebreak.m
   ================================================================
   peeling4Df4.m（Max Det タイブレーク版）を全6胞体・全 C1 に適用．
   旧版結果 ans4DGlobal.mx と比較し，ans4DGlobal_v2.mx に保存する．
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];
Get[scriptDir <> "peeling4Df4.m"];

dataDir = scriptDir <> "_4DData/";

polytopeNames = <|
  5   -> "5-cell",
  8   -> "8-cell",
  16  -> "16-cell",
  24  -> "24-cell",
  120 -> "120-cell",
  600 -> "600-cell"
|>;

loadPoly[n_] := Module[{raw},
  raw = Get[dataDir <> "f" <> ToString[n] <> ".m"];
  <|"vers" -> raw[[1]], "faces" -> raw[[4]], "cells" -> raw[[6]]|>
];

(* 旧版結果を読み込む *)
oldFile = scriptDir <> "ans4DGlobal.mx";
If[FileExistsQ[oldFile], Get[oldFile]; hasOld = True,
   Print["WARNING: ans4DGlobal.mx not found."]; hasOld = False];

Print[""];
Print[StringRepeat["=", 65]];
Print["  Apple-Peel 4D: all C1, all C2 (Det tie-breaking)"];
Print[StringRepeat["=", 65]];

results4DGlobal_v2 = <||>;

Do[
  poly = loadPoly[n];
  vers = poly["vers"]; faces = poly["faces"]; cells = poly["cells"];
  nc   = Length[cells];
  Print["\nProcessing ", polytopeNames[n], " (", nc, " cells) ..."];
  t0 = AbsoluteTime[];

  ans = Table[peeling4Df4[N[vers], faces, cells, top], {top, nc}];

  dt         = Round[AbsoluteTime[] - t0, 0.1];
  allBools   = Flatten[Map[Last, ans, {2}]];
  trueCount  = Count[allBools, True];
  falseCount = Count[allBools, False];
  uniqOrders = DeleteDuplicates[
    Select[Flatten[ans, 1], Last[#]&][[All, 2]]];
  classif = Which[falseCount == 0, "Perfect",
                  trueCount  == 0, "Impossible",
                  True,            "Possible"];

  Print["  Done: ", dt, " s  |  ", classif,
        "  (True: ", trueCount, ", False: ", falseCount,
        ", Unique: ", Length[uniqOrders], ")"];

  (* 旧版との比較 *)
  If[hasOld,
    oldTrue  = results4DGlobal[n]["trueCount"];
    oldFalse = results4DGlobal[n]["falseCount"];
    oldUniq  = results4DGlobal[n]["uniqueCount"];
    oldClass = results4DGlobal[n]["classification"];
    Print["  Old:  ", oldClass,
          "  (True: ", oldTrue, ", False: ", oldFalse,
          ", Unique: ", oldUniq, ")"];
    Print["  Diff: True ", oldTrue, " -> ", trueCount,
          If[trueCount > oldTrue, " (+", If[trueCount < oldTrue, " (", " (="]],
          If[trueCount != oldTrue, ToString[trueCount - oldTrue] <> ")", ")"],
          "  Unique ", oldUniq, " -> ", Length[uniqOrders]]
  ];

  results4DGlobal_v2[n] = <|
    "name"           -> polytopeNames[n],
    "n"              -> n,
    "classification" -> classif,
    "trueCount"      -> trueCount,
    "falseCount"     -> falseCount,
    "uniqueCount"    -> Length[uniqOrders],
    "results"        -> ans
  |>;
, {n, {5, 8, 16, 24, 120, 600}}];

(* 保存 *)
outFile = scriptDir <> "ans4DGlobal_v2.mx";
DumpSave[outFile, results4DGlobal_v2];
Print["\nSaved: ", outFile];

(* サマリ表 *)
Print[""];
Print[StringRepeat["=", 75]];
Print["  Summary (new vs old)"];
Print[StringRepeat["=", 75]];
Print[StringPadLeft["Polytope", 10], "  ",
      StringPadLeft["New Class", 10], "  ",
      StringPadLeft["True(n)", 8], "  ",
      StringPadLeft["True(o)", 8], "  ",
      StringPadLeft["Uniq(n)", 8], "  ",
      StringPadLeft["Uniq(o)", 8]];
Print[StringRepeat["-", 65]];
Do[
  r = results4DGlobal_v2[n];
  If[hasOld,
    o = results4DGlobal[n];
    Print[StringPadLeft[polytopeNames[n], 10], "  ",
          StringPadLeft[r["classification"], 10], "  ",
          StringPadLeft[ToString[r["trueCount"]],  8], "  ",
          StringPadLeft[ToString[o["trueCount"]],  8], "  ",
          StringPadLeft[ToString[r["uniqueCount"]], 8], "  ",
          StringPadLeft[ToString[o["uniqueCount"]], 8]],
    Print[StringPadLeft[polytopeNames[n], 10], "  ",
          StringPadLeft[r["classification"], 10], "  ",
          StringPadLeft[ToString[r["trueCount"]], 8]]
  ],
  {n, {5, 8, 16, 24, 120, 600}}
];

Print["\nDone."];
