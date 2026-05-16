(* ================================================================
   run_dodecahedron_epsilon.m
   正十二面体について左半空間条件の閾値 epsilon を変えて比較する．

   厳格 (strict):  Det[...] >= +epsilon  (境界近傍を左候補から除外)
   現行 (current): Det[...] >= 0
   寛容 (lenient): Det[...] >= -epsilon  (境界近傍も左候補に含める)

   各条件で R1/R2/R3 × fallback あり/なし の成功数を集計する．
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "peeling3DLoxo.m"];

(* --- 閾値を引数にとる汎用ペーリング関数 --- *)
p3lPeelPairEps[cc_, adj_, nf_, top_, f2_, rule_, eps_, useFallback_] :=
  Module[{order, inOrder, nbrs, prev, last,
          pending, leftCands, nHat, fHat, lHat,
          angles, absAngles, next},

    order   = {top, f2};
    inOrder = <|top -> True, f2 -> True|>;
    nbrs    = Table[DeleteCases[DeleteCases[adj[[i]], top], f2], {i, nf}];
    prev = top; last = f2;

    While[Length[order] < nf,
      pending = Select[nbrs[[last]], !KeyExistsQ[inOrder, #]&];
      nbrs[[last]] = pending;
      If[Length[pending] == 0, Break[]];

      {nHat, fHat, lHat} = p3lDarbouxFrame[cc[[prev]], cc[[last]]];

      (* 左半空間条件（閾値 eps を使用）
         候補が1つだけの場合は左条件によらず選択する *)
      leftCands = If[Length[pending] == 1,
        pending,
        Select[pending, Det[{cc[[prev]], cc[[last]], cc[[#]]}] >= eps &]
      ];

      next = Which[
        rule === "maxphi" && Length[leftCands] > 0,
          angles = p3lAngleFromForward[fHat, lHat, cc[[#]]] & /@ leftCands;
          leftCands[[ First[Ordering[-angles, 1]] ]],

        rule === "loxo" && Length[leftCands] > 0,
          angles = p3lAngleFromForward[fHat, lHat, cc[[#]]] & /@ leftCands;
          leftCands[[ First[Ordering[angles,  1]] ]],

        rule === "maxz" && Length[leftCands] > 0,
          leftCands[[ First[Ordering[-cc[[leftCands, 3]], 1]] ]],

        (* フォールバック *)
        useFallback && rule === "maxphi" && Length[pending] > 0,
          absAngles = Abs[p3lAngleFromForward[fHat, lHat, cc[[#]]]] & /@ pending;
          pending[[ First[Ordering[-absAngles, 1]] ]],

        useFallback && rule === "loxo" && Length[pending] > 0,
          absAngles = Abs[p3lAngleFromForward[fHat, lHat, cc[[#]]]] & /@ pending;
          pending[[ First[Ordering[absAngles, 1]] ]],

        useFallback && rule === "maxz" && Length[pending] > 0,
          pending[[ First[Ordering[cc[[pending, 3]], 1]] ]],

        True, Break[]; Null
      ];

      If[next === Null, Break[]];
      AppendTo[order, next];
      inOrder[next] = True;
      nbrs = Table[DeleteCases[nbrs[[i]], next], {i, nf}];
      prev = last; last = next;
    ];
    {order, Length[order] == nf}
  ];

p3lRunAllEps[vers_, faces_, rule_, eps_, useFallback_] :=
  Module[{nf, adj, lv, cc, results, f2Cands},
    nf  = Length[faces];
    adj = p3lBuildAdj[faces];
    results = {};
    Do[
      lv      = p3lAlignTopToZ[vers, faces, top];
      cc      = p3lFaceCentroid[#, lv, faces] & /@ Range[nf];
      f2Cands = adj[[top]];
      Do[
        Module[{res},
          res = p3lPeelPairEps[cc, adj, nf, top, f2Cands[[fi]], rule, eps, useFallback];
          AppendTo[results, {top, f2Cands[[fi]], res[[1]], res[[2]]}]
        ],
        {fi, Length[f2Cands]}
      ],
      {top, nf}
    ];
    results
  ];

(* --- Dodecahedron データ読み込み --- *)
{vers, faces} = p3lExtractPolyhedron["Dodecahedron"];
nf  = Length[faces];
adj = p3lBuildAdj[faces];
totalPairs = Total[Length /@ adj];

Print[""];
Print[StringRepeat["=", 72]];
Print["  Dodecahedron: 左半空間条件の閾値比較 (総ペア数 = ", totalPairs, ")"];
Print[StringRepeat["=", 72]];

(* --- 閾値の一覧 --- *)
epsValues = {1*^-6, 1*^-8, 1*^-10, 1*^-12, 0, -1*^-12, -1*^-10, -1*^-8, -1*^-6};
epsLabels = {"strict 1e-6", "strict 1e-8", "strict 1e-10", "strict 1e-12",
             "current 0",
             "lenient 1e-12", "lenient 1e-10", "lenient 1e-8", "lenient 1e-6"};

rules      = {"maxphi", "loxo", "maxz"};
ruleLabels = {"R1", "R2", "R3"};

(* ヘッダ *)
Print[""];
Print[StringPadLeft["epsilon", 16], "   ",
      StringPadLeft["R1(w)", 7],
      StringPadLeft["R1(n)", 7],
      StringPadLeft["R2(w)", 7],
      StringPadLeft["R2(n)", 7],
      StringPadLeft["R3(w)", 7],
      StringPadLeft["R3(n)", 7]];
Print[StringRepeat["-", 70]];

Do[
  eps   = epsValues[[ei]];
  label = epsLabels[[ei]];
  row   = label;
  cols  = {};

  Do[
    rw = p3lRunAllEps[vers, faces, rules[[ri]], eps, True];
    nw = Length[Select[rw, #[[4]]&]];
    rn = p3lRunAllEps[vers, faces, rules[[ri]], eps, False];
    nn = Length[Select[rn, #[[4]]&]];
    AppendTo[cols, {nw, nn}],
    {ri, 3}
  ];

  Print[StringPadLeft[label, 16], "   ",
        Sequence @@ Table[
          StringPadLeft[ToString[cols[[ri, 1]]], 7] <>
          StringPadLeft[ToString[cols[[ri, 2]]], 7],
          {ri, 3}
        ]
  ],
  {ei, Length[epsValues]}
];

Print[StringRepeat["-", 70]];
Print["(w) = with fallback,  (n) = no fallback"];
Print["単一候補優先は全ケースに適用"];
Print["\nDone."];
