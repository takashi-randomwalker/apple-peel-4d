(* ================================================================
   compare_TI_R1_R3.m
   切頂二十面体における R1 (maxphi) vs R3 (maxz) の比較
   peeling3DLoxo.m の基盤を流用し R1 を追加
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "peeling3DLoxo.m"];

(* --- R1 "maxphi" 規則を追加した 1ペア実行関数 --- *)
(* p3lPeelPair は "loxo" / "maxz" のみ対応なので
   3規則対応の wrapper を定義する *)
p3lPeelPairR1[cc_, adj_, nf_, top_, f2_] :=
  Module[{order, inOrder, nbrs, prev, last,
          pending, lcands, frame, nHat, fHat, lHat,
          angles, absAngles, idx, next},

    order   = {top, f2};
    inOrder = <|top -> True, f2 -> True|>;
    nbrs    = Table[
      DeleteCases[DeleteCases[adj[[i]], top], f2],
      {i, nf}
    ];
    prev = top;
    last = f2;

    While[Length[order] < nf,
      pending = Select[nbrs[[last]], !KeyExistsQ[inOrder, #]&];
      nbrs[[last]] = pending;
      If[Length[pending] == 0, Break[]];

      {nHat, fHat, lHat} = p3lDarbouxFrame[cc[[prev]], cc[[last]]];
      lcands = Select[pending,
        Det[{cc[[prev]], cc[[last]], cc[[#]]}] >= 0 &];

      If[Length[lcands] > 0,
        (* R1: 左候補から max φ *)
        angles = p3lAngleFromForward[fHat, lHat, cc[[#]]] & /@ lcands;
        next = lcands[[ First[Ordering[-angles, 1]] ]],
        (* フォールバック: 全候補から max |φ| *)
        absAngles = Abs[p3lAngleFromForward[fHat, lHat, cc[[#]]]] & /@ pending;
        next = pending[[ First[Ordering[-absAngles, 1]] ]]
      ];

      AppendTo[order, next];
      inOrder[next] = True;
      nbrs = Map[DeleteCases[#, next]&, nbrs];
      prev = last; last = next;
    ];
    {order, Length[order] == nf}
  ];

(* 全ペア実行 (R1) *)
p3lRunAllR1[vers_, faces_] :=
  Module[{nf, adj, lv, cc, results, f2Cands, order, complete},
    nf      = Length[faces];
    adj     = p3lBuildAdj[faces];
    results = {};
    Do[
      lv  = p3lAlignTopToZ[vers, faces, top];
      cc  = p3lFaceCentroid[#, lv, faces] & /@ Range[nf];
      f2Cands = adj[[top]];
      Do[
        {order, complete} = p3lPeelPairR1[cc, adj, nf, top, f2Cands[[fi]]];
        AppendTo[results, {top, f2Cands[[fi]], order, complete}],
        {fi, Length[f2Cands]}
      ],
      {top, nf}
    ];
    results
  ];

(* ----------------------------------------------------------------
   切頂二十面体データの読み込み
   ---------------------------------------------------------------- *)
Print["Loading TruncatedIcosahedron..."];
{vers, faces} = p3lExtractPolyhedron["TruncatedIcosahedron"];
nf = Length[faces];
adj = p3lBuildAdj[faces];
Print["Faces: ", nf, "  Total pairs: ", Total[Length /@ adj]];

(* 面種類: 頂点数 5 = 五角形, 6 = 六角形 *)
faceType[fi_] := If[Length[faces[[fi]]] == 5, "P", "H"];

(* ----------------------------------------------------------------
   R1 / R3 を全ペア実行
   ---------------------------------------------------------------- *)
Print["Running R1 (maxphi)..."];
r1All = p3lRunAllR1[vers, faces];
r1Suc = Select[r1All, #[[4]]&];

Print["Running R3 (maxz)..."];
r3All = p3lRunAll[vers, faces, "maxz"];
r3Suc = Select[r3All, #[[4]]&];

Print[""];
Print["=== 切頂二十面体 R1 vs R3 ==="];
Print["R1 success: ", Length[r1Suc], " / ", Length[r1All],
      "  (", N[100 Length[r1Suc]/Length[r1All], 4], "%)"];
Print["R3 success: ", Length[r3Suc], " / ", Length[r3All],
      "  (", N[100 Length[r3Suc]/Length[r3All], 4], "%)"];

(* ----------------------------------------------------------------
   成功集合の比較
   ---------------------------------------------------------------- *)
r1Keys = {#[[1]], #[[2]]} & /@ r1Suc;
r3Keys = {#[[1]], #[[2]]} & /@ r3Suc;

r1Only   = Select[r1Suc, !MemberQ[r3Keys, {#[[1]],#[[2]]}]&];
r3Only   = Select[r3Suc, !MemberQ[r1Keys, {#[[1]],#[[2]]}]&];
bothSucc = Select[r1Suc,  MemberQ[r3Keys, {#[[1]],#[[2]]}]&];

Print[""];
Print["R1のみ成功: ", Length[r1Only]];
Print["R3のみ成功: ", Length[r3Only]];
Print["両方成功:   ", Length[bothSucc]];

(* ----------------------------------------------------------------
   代表ペアの face sequence & 面種別パターン
   ---------------------------------------------------------------- *)
showR1only = Take[r1Only,   Min[3, Length[r1Only]]];
showBoth   = Take[bothSucc, Min[2, Length[bothSucc]]];
showR3only = Take[r3Only,   Min[2, Length[r3Only]]];

printPair[top_, f2_] :=
  Module[{r1m, r3m, seqR1, seqR3, patR1, patR3},
    r1m = SelectFirst[r1All, #[[1]] == top && #[[2]] == f2 &];
    r3m = SelectFirst[r3All, #[[1]] == top && #[[2]] == f2 &];
    seqR1 = If[MissingQ[r1m], {}, r1m[[3]]];
    seqR3 = If[MissingQ[r3m], {}, r3m[[3]]];
    patR1 = StringJoin[faceType /@ seqR1];
    patR3 = StringJoin[faceType /@ seqR3];
    Print["  (F1,F2)=(", top, ",", f2, ")"];
    Print["  R1 seq : ", seqR1, "  [", Length[seqR1], "/", nf, "  ", If[Length[seqR1]==nf,"OK","NG"], "]"];
    Print["  R3 seq : ", seqR3, "  [", Length[seqR3], "/", nf, "  ", If[Length[seqR3]==nf,"OK","NG"], "]"];
    Print["  R1 pat : ", patR1];
    Print["  R3 pat : ", patR3];
    Print[""];
  ];

Print[""];
Print["--- R1 のみ成功のペア ---"];
printPair[#[[1]], #[[2]]] & /@ showR1only;

Print["--- 両方成功のペア ---"];
printPair[#[[1]], #[[2]]] & /@ showBoth;

Print["--- R3 のみ成功のペア ---"];
printPair[#[[1]], #[[2]]] & /@ showR3only;

(* ----------------------------------------------------------------
   面種別パターン統計: P/H の交互パターン分析
   ---------------------------------------------------------------- *)
(* 連続する同種面の最大ラン長 = ジグザグ度の指標 *)
maxRun[pat_] :=
  Max[Length /@ SplitBy[Characters[pat], Identity]];

If[Length[r1Suc] > 0,
  r1MaxRuns = maxRun[StringJoin[faceType /@ #[[3]]]] & /@ r1Suc;
  Print["R1成功例の max run length (小ほど交互): mean=",
        N[Mean[r1MaxRuns], 3], " range=", {Min[r1MaxRuns], Max[r1MaxRuns]}];
];
If[Length[r3Suc] > 0,
  r3MaxRuns = maxRun[StringJoin[faceType /@ #[[3]]]] & /@ r3Suc;
  Print["R3成功例の max run length (小ほど交互): mean=",
        N[Mean[r3MaxRuns], 3], " range=", {Min[r3MaxRuns], Max[r3MaxRuns]}];
];

Print["Done."];
