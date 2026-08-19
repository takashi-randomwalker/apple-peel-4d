(* ================================================================
   run4DBacktrack.m
   ================================================================
   peeling4Df4_1back.m（1段バックトラック，RZ 規則）を
   全6正多胞体・全 C1 に適用し，グリーディ版（ans4DGlobal_v3.mx）と比較する．
   結果を ans4DGlobal1back.mx に保存する．
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];
Get[scriptDir <> "peeling4Df4_1back.m"];

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

summarize[ans_] := Module[
  {allBools, trueCount, falseCount, uniqOrders, classif},
  allBools   = Flatten[Map[Last, ans, {2}]];
  trueCount  = Count[allBools, True];
  falseCount = Count[allBools, False];
  uniqOrders = DeleteDuplicates[
    Select[Flatten[ans, 1], Last[#] &][[All, 2]]
  ];
  classif = Which[
    falseCount == 0, "Perfect",
    trueCount  == 0, "Impossible",
    True,            "Possible"
  ];
  <|"trueCount"      -> trueCount,
    "falseCount"     -> falseCount,
    "uniqueCount"    -> Length[uniqOrders],
    "classification" -> classif|>
];

(* グリーディ版結果（比較用）：DownValues 形式で格納されている *)
v3File = scriptDir <> "ans4DGlobal_v3.mx";
If[FileExistsQ[v3File],
  Get[v3File]; hasV3 = True;
  Print["Loaded greedy results from ans4DGlobal_v3.mx"],
  Print["WARNING: ans4DGlobal_v3.mx not found. Comparison skipped."];
  hasV3 = False
];

Print[""];
Print[StringRepeat["=", 65]];
Print["  Apple-Peel 4D: 1-step backtrack (RZ), all C1/C2"];
Print[StringRepeat["=", 65]];

(* 結果をシンボルの DownValues として格納（下線なし） *)
Clear[results4D1back];

Do[
  poly  = loadPoly[n];
  vers  = poly["vers"];
  faces = poly["faces"];
  cells = poly["cells"];
  nc    = Length[cells];
  Print["\n--- ", polytopeNames[n], " (", nc, " cells) ---"];
  t0 = AbsoluteTime[];

  ans = Table[
    If[Mod[top, 100] == 0 && nc >= 100,
      Print["  C1=", top, " / ", nc,
            "  (", Round[AbsoluteTime[] - t0, 1], " s elapsed)"]
    ];
    peeling4Df41back[N[vers], faces, cells, top, "RZ"],
    {top, nc}
  ];

  dt = Round[AbsoluteTime[] - t0, 0.1];
  s  = summarize[ans];

  Print["  1back(RZ): ", dt, " s  |  ", s["classification"],
        "  (True: ", s["trueCount"],
        ", False: ", s["falseCount"],
        ", Unique: ", s["uniqueCount"], ")"];

  (* ValueQ でダウンバリューの存在を確認 *)
  If[hasV3 && ValueQ[results4DGlobal_v3[n]],
    v    = results4DGlobal_v3[n];
    diff = s["trueCount"] - v["trueCount"];
    Print["  greedy(RZ): ", v["classification"],
          "  (True: ", v["trueCount"],
          ", Unique: ", v["uniqueCount"], ")"];
    Print["  Improvement: ",
          If[diff > 0, "+", ""], diff, " True pairs"]
  ];

  results4D1back[n] = Join[s, <|
    "name"    -> polytopeNames[n],
    "n"       -> n,
    "results" -> ans
  |>],
  {n, {5, 8, 16, 24, 120, 600}}
];

(* 保存 *)
outFile = scriptDir <> "ans4DGlobal1back.mx";
DumpSave[outFile, results4D1back];
Print["\nSaved: ", outFile];

(* ── サマリ表 ─────────────────────────────────────────────── *)
Print[""];
Print[StringRepeat["=", 75]];
Print["  Summary: 1-step backtrack (RZ) vs greedy (RZ)  [ans4DGlobal_v3]"];
Print[StringRepeat["=", 75]];
Print[StringPadLeft["Polytope",  10], "  ",
      StringPadLeft["1back",     10], "  ",
      StringPadLeft["True(1b)",   9], "  ",
      StringPadLeft["True(gr)",   9], "  ",
      StringPadLeft["Diff",       6], "  ",
      StringPadLeft["Uniq(1b)",   9]];
Print[StringRepeat["-", 62]];
Do[
  r = results4D1back[n];
  If[hasV3 && ValueQ[results4DGlobal_v3[n]],
    v    = results4DGlobal_v3[n];
    diff = r["trueCount"] - v["trueCount"];
    Print[StringPadLeft[polytopeNames[n],        10], "  ",
          StringPadLeft[r["classification"],      10], "  ",
          StringPadLeft[ToString[r["trueCount"]],  9], "  ",
          StringPadLeft[ToString[v["trueCount"]],  9], "  ",
          StringPadLeft[If[diff >= 0,
                           "+" <> ToString[diff],
                           ToString[diff]],        6], "  ",
          StringPadLeft[ToString[r["uniqueCount"]], 9]],
    Print[StringPadLeft[polytopeNames[n],          10], "  ",
          StringPadLeft[r["classification"],        10], "  ",
          StringPadLeft[ToString[r["trueCount"]],    9]]
  ],
  {n, {5, 8, 16, 24, 120, 600}}
];

Print["\nDone."];
