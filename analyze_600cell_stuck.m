(* ================================================================
   analyze_600cell_stuck.m
   ================================================================
   600 胞体の全 (C1,C2) ペアを peeling4Df4 で実行し，
   終了時の order 長（= 配置済みセル数）の頻度分布を表示する．
   localVers は保存せず，order 長のみ収集する．
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];
Get[scriptDir <> "peeling4Df4.m"];

dataDir = scriptDir <> "_4DData/";
raw   = Get[dataDir <> "f600.m"];
vers  = raw[[1]];
faces = raw[[4]];
cells = raw[[6]];
nc    = Length[cells];  (* 600 *)

Print["600-cell: ", nc, " cells, computing all (C1,C2) pairs..."];
t0 = AbsoluteTime[];

(* 全 C1 について peeling を実行し，終了ステップ長のみ収集 *)
stuckLengths = {};
Do[
  res = peeling4Df4[N[vers], faces, cells, top];
  (* res = { {localVers, order, success}, ... } ← C2 候補数分 *)
  lengths = Length[#[[2]]] & /@ res;
  stuckLengths = Join[stuckLengths, lengths];
  If[Mod[top, 100] == 0,
     Print["  C1=", top, " done (", Length[stuckLengths], " pairs so far, ",
           Round[AbsoluteTime[] - t0, 1], " s)"]],
  {top, nc}
];

dt = Round[AbsoluteTime[] - t0, 0.1];
Print["Computation done: ", dt, " s"];
Print["Total pairs: ", Length[stuckLengths]];
Print[""];

successes = Count[stuckLengths, nc];
failures  = Length[stuckLengths] - successes;
Print["Successes (order length = ", nc, "): ", successes];
Print["Failures:  ", failures];
Print[""];

(* 頻度集計 *)
tally = Sort[Tally[stuckLengths]];
maxCount = Max[tally[[All, 2]]];
total    = Length[stuckLengths];

Print[StringRepeat["=", 58]];
Print["  Termination step distribution (600-cell, ", total, " pairs)"];
Print[StringRepeat["=", 58]];
Print[StringPadLeft["Step", 6], "  ",
      StringPadLeft["Count", 6], "  ",
      StringPadLeft["%", 6], "  Bar (max=", maxCount, ")"];
Print[StringRepeat["-", 58]];
Do[
  step = t[[1]];
  cnt  = t[[2]];
  pct  = N[cnt / total * 100, 3];
  barW = Max[1, Round[cnt / maxCount * 40]];
  Print[
    StringPadLeft[ToString[step], 6], "  ",
    StringPadLeft[ToString[cnt], 6], "  ",
    StringPadLeft[ToString[PaddedForm[pct, {4, 1}]], 6], "  ",
    StringRepeat["#", barW]
  ],
  {t, tally}
];
Print[StringRepeat["-", 58]];

(* 統計サマリ *)
Print[""];
Print["Min stuck step: ", Min[stuckLengths]];
Print["Max stuck step: ", Max[stuckLengths]];
Print["Mean:           ", N[Mean[stuckLengths], 4]];
Print["Median:         ", Median[stuckLengths]];

(* 累積分布（代表値のみ） *)
Print[""];
Print[StringRepeat["=", 45]];
Print["  Cumulative distribution"];
Print[StringRepeat["=", 45]];
cumSum = 0;
prevStep = 0;
Do[
  step = t[[1]];
  cnt  = t[[2]];
  cumSum += cnt;
  cumPct = N[cumSum / total * 100, 3];
  Print["  <= ", StringPadLeft[ToString[step], 4], ": ",
        StringPadLeft[ToString[PaddedForm[cumPct, {5, 1}]], 6],
        "%  (", cumSum, ")"],
  {t, tally}
];

Print["\nDone."];
