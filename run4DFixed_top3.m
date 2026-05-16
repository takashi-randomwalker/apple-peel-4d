(* ================================================================
   run4DFixed_top3.m
   ================================================================
   peeling4Df4.m（タイブレーク改訂版）を C1 = 3 固定で実行し，
   全6胞体について隣接するすべての C2 を調べる．

   比較: ans4DGlobal.mx（タイブレーク前の旧版計算結果）と
         新版（Max Det タイブレーク）の成否の差異を報告する．
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];
Get[scriptDir <> "peeling4Df4.m"];

dataDir = scriptDir <> "_4DData/";

fixedTop = 3;

polytopeNames = <|
  5   -> "5-cell",
  8   -> "8-cell",
  16  -> "16-cell",
  24  -> "24-cell",
  120 -> "120-cell",
  600 -> "600-cell"
|>;

loadPoly[n_] := Module[{raw, vers, faces, cells},
  raw   = Get[dataDir <> "f" <> ToString[n] <> ".m"];
  vers  = raw[[1]];
  faces = raw[[4]];
  cells = raw[[6]];
  <|"vers" -> vers, "faces" -> faces, "cells" -> cells|>
];

(* 旧版結果を読み込む *)
oldFile = scriptDir <> "ans4DGlobal.mx";
hasOld = FileExistsQ[oldFile];
If[hasOld,
  Get[oldFile];
  Print["Loaded old results from ans4DGlobal.mx"],
  Print["WARNING: ans4DGlobal.mx not found; skipping comparison."]
];

Print[""];
Print[StringRepeat["=", 65]];
Print["  Apple-Peel 4D: fixed C1=", fixedTop, ", all adjacent C2"];
Print["  Algorithm: peeling4Df4 + Det tie-breaking"];
Print[StringRepeat["=", 65]];

newResults = <||>;

Do[
  poly = loadPoly[n];
  vers  = poly["vers"];
  faces = poly["faces"];
  cells = poly["cells"];
  nc    = Length[cells];

  t0  = AbsoluteTime[];
  res = peeling4Df4[N[vers], faces, cells, fixedTop];
  dt  = AbsoluteTime[] - t0;

  (* 新版の成否リスト *)
  newBools  = Last /@ res;
  nC2total  = Length[res];
  nSuccess  = Count[newBools, True];

  (* 旧版の成否リスト（C1=fixedTop のみ取り出す） *)
  If[hasOld,
    oldRes   = results4DGlobal[n]["results"][[fixedTop]];
    oldBools = Last /@ oldRes;
    same = (newBools === oldBools);
    nChanged = Count[Thread[newBools != oldBools], True],
    same = True; nChanged = 0
  ];

  (* 隣接 C2 インデックス *)
  adj = Select[Range[nc], Function[j,
    j != fixedTop &&
    Length[Intersection[cells[[fixedTop]], cells[[j]]]] >= 1]];

  Print["\n--- ", polytopeNames[n],
        "  (cells=", nc, ", C1=", fixedTop, ", #C2=", nC2total, ")",
        "  ", NumberForm[dt, {4,2}], " s ---"];
  Do[
    c2idx = adj[[i]];
    nb = newBools[[i]];
    If[hasOld,
      ob = oldBools[[i]];
      tag = If[nb === ob, "  ", If[nb && !ob, " NEW-OK", " NEW-FAIL"]];
      Print["  C2=", c2idx, "  new=", nb, "  old=", ob, tag],
      Print["  C2=", c2idx, "  success=", nb]
    ],
    {i, nC2total}
  ];

  Print["  >> new: ", nSuccess, "/", nC2total, " success"];
  If[hasOld,
    Print["  >> old: ", Count[oldBools, True], "/", nC2total, " success"];
    Print["  >> same: ", same, "  changed: ", nChanged, " pair(s)"]
  ];

  newResults[n] = <|"bools" -> newBools, "adj" -> adj, "res" -> res|>;

, {n, {5, 8, 16, 24, 120, 600}}];

Print[""];
Print[StringRepeat["=", 65]];
Print["  Summary (C1=", fixedTop, ")"];
Print[StringRepeat["=", 65]];
Print[StringPadLeft["Polytope", 10], "  ",
      StringPadLeft["C2s", 4], "  ",
      StringPadLeft["New", 4], "  ",
      If[hasOld, StringPadLeft["Old", 4] <> "  " <> StringPadLeft["Diff", 4], ""]];
Print[StringRepeat["-", If[hasOld, 35, 25]]];
Do[
  r   = newResults[n];
  nb  = r["bools"];
  tot = Length[nb];
  ns  = Count[nb, True];
  If[hasOld,
    ob  = Last /@ results4DGlobal[n]["results"][[fixedTop]];
    os  = Count[ob, True];
    nd  = Count[Thread[nb != ob], True];
    Print[StringPadLeft[polytopeNames[n], 10], "  ",
          StringPadLeft[ToString[tot], 4], "  ",
          StringPadLeft[ToString[ns], 4], "  ",
          StringPadLeft[ToString[os], 4], "  ",
          StringPadLeft[ToString[nd], 4]],
    Print[StringPadLeft[polytopeNames[n], 10], "  ",
          StringPadLeft[ToString[tot], 4], "  ",
          StringPadLeft[ToString[ns], 4]]
  ],
  {n, {5, 8, 16, 24, 120, 600}}
];

Print["\nDone."];
