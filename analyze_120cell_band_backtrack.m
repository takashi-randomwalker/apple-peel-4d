(* ================================================================
   analyze_120cell_band_backtrack.m
   帯境界の代替エントリーセル選択で自己交差を削減できるか検証
   C1=3, C2=1 (Type A 代表) を対象とする
   各帯遷移ステップで，RZ が選ぶセルの代わりに隣接する同帯セルを試し
   重なり数が減るか（0 になるか）を調べる
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName], DirectoryName[$InputFileName], NotebookDirectory[]];
Get[scriptDir <> "unfold3DExport.m"];
Get[scriptDir <> "ans4DGlobal_v4.mx"];

dataDir = scriptDir <> "_4DData/";
raw   = Get[dataDir <> "f120.m"];
vers  = N[raw[[1]]]; faces = raw[[4]]; cells = raw[[6]];
nc    = Length[cells];

cellVerts[ci_] := Union @@ faces[[ cells[[ci]] ]];
eps    = 10^-10;
detEps = 10^-10;

(* ----------------------------------------------------------------
   C1=3, C2=1 の結果取得
   ---------------------------------------------------------------- *)
fixedC1  = 3; targetC2 = 1;
allRes   = results4DGlobal[120]["results"][[fixedC1]];
rep      = SelectFirst[allRes, Last[#] && #[[2, 2]] == targetC2 &];
localVers = rep[[1]];
baseOrd   = rep[[2]];

(* ----------------------------------------------------------------
   セル重心とバンド割り当て
   ---------------------------------------------------------------- *)
cC      = Mean[localVers[[ cellVerts[#] ]]] & /@ Range[nc];
wRound  = Round[cC[[#, 4]] & /@ Range[nc], 0.01];
wLevels = Sort[DeleteDuplicates[wRound], Greater];
bandOf  = Association @ Table[c -> Position[wLevels, wRound[[c]]][[1, 1]], {c, nc}];

c1Vec = cC[[fixedC1]];
c2Vec = cC[[targetC2]];

(* ----------------------------------------------------------------
   全隣接リスト（一度だけ計算）
   ---------------------------------------------------------------- *)
Print["Building adjacency list..."];
t0 = AbsoluteTime[];
fullAdj = Table[
  Select[Range[nc], j |-> j != i &&
    Length[Intersection[cells[[i]], cells[[j]]]] >= 1],
  {i, nc}];
Print["  Done: ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* ----------------------------------------------------------------
   RZ グリーディ継続関数
   partial: 確定済み部分順序（先頭が C1，2番目が C2）
   戻り値: {完全順序, 成功True/False}
   ---------------------------------------------------------------- *)
continueRZ[partial_] :=
  Module[{order, adjList, last, cans, rCans, wvals, maxw, tied,
          gvals, xyvals, newCell, gs, xs},
    order   = partial;
    adjList = Complement[fullAdj[[#]], order] & /@ Range[nc];
    last    = Last[order];
    gs[k_] := If[Length[order] == 2,
      c2Vec[[1]] * cC[[k, 2]] - c2Vec[[2]] * cC[[k, 1]],
      Quiet[Det[{c1Vec, c2Vec, cC[[last]], cC[[k]]}]]];
    xs[k_] := c2Vec[[1]] * cC[[k, 2]] - c2Vec[[2]] * cC[[k, 1]];
    While[Length[Flatten[adjList]] != 0 && Length[adjList[[last]]] != 0,
      cans = If[Length[order] == 2, adjList[[last]],
        Select[adjList[[last]], j |->
          Quiet[Det[{c1Vec, c2Vec, cC[[last]], cC[[j]]}]] >= -eps]];
      newCell = If[Length[cans] >= 1,
        wvals = cC[[cans, 4]]; maxw = Max[wvals];
        tied = Select[cans, cC[[#, 4]] >= maxw - 10^-10 &];
        If[Length[tied] == 1, tied[[1]],
          gvals = gs /@ tied;
          If[Length[order] >= 3 && Max[Abs[gvals]] < detEps,
            xyvals = xs /@ tied; tied[[First[Ordering[-xyvals, 1]]]],
            tied[[First[Ordering[-gvals, 1]]]]]],
        rCans = adjList[[last]];
        wvals = cC[[rCans, 4]]; maxw = Min[wvals];
        tied = Select[rCans, cC[[#, 4]] <= maxw + 10^-10 &];
        If[Length[tied] == 1, tied[[1]],
          gvals = Abs[gs[#]] & /@ tied;
          If[Length[order] >= 3 && Max[gvals] < detEps,
            xyvals = Abs[xs[#]] & /@ tied; tied[[First[Ordering[-xyvals, 1]]]],
            tied[[First[Ordering[-gvals, 1]]]]]]];
      last = newCell;
      AppendTo[order, newCell];
      adjList = Complement[#, {newCell}] & /@ adjList;
    ];
    {order, Length[order] == nc}
  ];

(* ----------------------------------------------------------------
   cellLocal3Ds 事前計算
   ---------------------------------------------------------------- *)
Print["Pre-computing cellLocal3Ds..."];
t0 = AbsoluteTime[];
cellLocal3Ds = Table[unfoldCellLocal3D[ci, vers, cells, faces], {ci, nc}];
Print["  Done: ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* ----------------------------------------------------------------
   ベースライン重なり
   ---------------------------------------------------------------- *)
baseEmbs = unfoldTo3DFast[cells, faces, baseOrd, cellLocal3Ds];
baseOv   = Length[checkOverlaps[baseOrd, cells, faces, baseEmbs]];
Print["Baseline overlaps (C1=3, C2=1): ", baseOv];

(* ----------------------------------------------------------------
   帯遷移ステップ特定
   ---------------------------------------------------------------- *)
ordBands   = bandOf /@ baseOrd;
transSteps = Select[Range[2, nc], ordBands[[#]] != ordBands[[# - 1]] &];
Print["帯遷移ステップ: ", transSteps,
      "  帯: ", ordBands[[transSteps - 1]], " → ", ordBands[[transSteps]]];

(* ----------------------------------------------------------------
   各帯遷移で代替エントリーセルを試す
   ---------------------------------------------------------------- *)
Print["\n== 各帯遷移での代替エントリー検証 =="];
bestOv  = baseOv;
bestOrd = baseOrd;
bestDesc = "baseline";

Do[
  t        = transSteps[[bi]];
  prevB    = ordBands[[t - 1]];
  newB     = ordBands[[t]];
  partial  = baseOrd[[1 ;; t - 1]];
  lastCell = Last[partial];

  (* newB に属し，lastCell に隣接するセル（baseOrd[[t]] を含む全候補） *)
  alts = Select[fullAdj[[lastCell]],
    c |-> bandOf[c] == newB && !MemberQ[partial, c] &];

  If[Length[alts] <= 1,
    Print["帯", prevB, "→帯", newB, "  (step ", t, ")  代替なし"],
    Print["帯", prevB, "→帯", newB, "  (step ", t, ")",
          "  baseline=", baseOrd[[t]],
          "  代替候補=", Complement[alts, {baseOrd[[t]]}]];
    Do[
      altCell = alts[[a]];
      If[altCell == baseOrd[[t]], Continue[]];
      t1 = AbsoluteTime[];
      {newOrd, ok} = continueRZ[Append[partial, altCell]];
      If[!ok,
        Print["  ", altCell, " → FAIL (incomplete)"],
        embs = unfoldTo3DFast[cells, faces, newOrd, cellLocal3Ds];
        ov   = Length[checkOverlaps[newOrd, cells, faces, embs]];
        dt   = Round[AbsoluteTime[] - t1, 0.1];
        star = Which[ov == 0, "  ★★ ZERO!", ov < baseOv, "  ★改善", True, ""];
        Print["  alt=", altCell, "  ov=", ov, star, "  (", dt, " s)"];
        If[ov < bestOv, bestOv = ov; bestOrd = newOrd;
           bestDesc = "帯" <> ToString[prevB] <> "→" <> ToString[newB] <>
                      " alt=" <> ToString[altCell]]
      ],
      {a, Length[alts]}
    ]
  ],
  {bi, Length[transSteps]}
];

(* ----------------------------------------------------------------
   サマリ
   ---------------------------------------------------------------- *)
Print["\n========================================"];
Print["Baseline: ", baseOv, " overlaps"];
Print["Best:     ", bestOv, " overlaps  (", bestDesc, ")"];
If[bestOv < baseOv,
  Print["改善: ", baseOv - bestOv, " ペア削減 (",
        Round[100. (baseOv - bestOv) / baseOv, 1], "%)"],
  Print["改善なし"]
];
If[bestOv == 0, Print["★★ 重なりゼロの net が見つかりました！"]];
