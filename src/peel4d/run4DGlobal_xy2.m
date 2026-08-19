(* ================================================================
   run4DGlobal_xy2.m
   ================================================================
   peeling4Df4.m（k=2 xy クロス積幾何タイブレーク版）を全6胞体・全 C1 に適用．
   前回結果 ans4DGlobal_v2.mx（Det タイブレーク・k=2 はリスト順）と比較し，
   ans4DGlobal_v3.mx に保存する．
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

(* v2 結果を読み込む（比較用） *)
oldFile = scriptDir <> "ans4DGlobal_v2.mx";
If[FileExistsQ[oldFile], Get[oldFile]; hasOld = True,
   Print["WARNING: ans4DGlobal_v2.mx not found."]; hasOld = False];

Print[""];
Print[StringRepeat["=", 65]];
Print["  Apple-Peel 4D: all C1, all C2 (k=2 xy cross-product tie-break)"];
Print[StringRepeat["=", 65]];

results4DGlobal_v3 = <||>;

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

  (* v2 との比較 *)
  If[hasOld,
    oldTrue  = results4DGlobal_v2[n]["trueCount"];
    oldFalse = results4DGlobal_v2[n]["falseCount"];
    oldUniq  = results4DGlobal_v2[n]["uniqueCount"];
    oldClass = results4DGlobal_v2[n]["classification"];
    Print["  v2:   ", oldClass,
          "  (True: ", oldTrue, ", False: ", oldFalse,
          ", Unique: ", oldUniq, ")"];
    Print["  Diff: True ", oldTrue, " -> ", trueCount,
          If[trueCount > oldTrue, " (+", If[trueCount < oldTrue, " (", " (="]],
          If[trueCount != oldTrue, ToString[trueCount - oldTrue] <> ")", ")"],
          "  Unique ", oldUniq, " -> ", Length[uniqOrders]]
  ];

  results4DGlobal_v3[n] = <|
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
outFile = scriptDir <> "ans4DGlobal_v3.mx";
DumpSave[outFile, results4DGlobal_v3];
Print["\nSaved: ", outFile];

(* サマリ表 *)
Print[""];
Print[StringRepeat["=", 75]];
Print["  Summary (v3: xy-cross k=2  vs  v2: list-order k=2)"];
Print[StringRepeat["=", 75]];
Print[StringPadLeft["Polytope", 10], "  ",
      StringPadLeft["Class(v3)", 10], "  ",
      StringPadLeft["True(v3)", 8], "  ",
      StringPadLeft["True(v2)", 8], "  ",
      StringPadLeft["Uniq(v3)", 8], "  ",
      StringPadLeft["Uniq(v2)", 8]];
Print[StringRepeat["-", 65]];
Do[
  r = results4DGlobal_v3[n];
  If[hasOld,
    o = results4DGlobal_v2[n];
    Print[StringPadLeft[r["name"], 10], "  ",
          StringPadLeft[r["classification"], 10], "  ",
          StringPadLeft[ToString[r["trueCount"]],   8], "  ",
          StringPadLeft[ToString[o["trueCount"]],   8], "  ",
          StringPadLeft[ToString[r["uniqueCount"]], 8], "  ",
          StringPadLeft[ToString[o["uniqueCount"]], 8]],
    Print[StringPadLeft[r["name"], 10], "  ",
          StringPadLeft[r["classification"], 10], "  ",
          StringPadLeft[ToString[r["trueCount"]], 8]]
  ],
  {n, {5, 8, 16, 24, 120, 600}}
];

Print["\nDone."];
