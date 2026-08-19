(* ================================================================
   run_platonic_eps_comparison.m
   左条件 det(c_1, c_k, c_j) >= thresh における閾値の影響を調べる

   3種類の閾値を比較:
     thresh = -epsVal  (eps > 0: 現在の設定, 許容あり)
     thresh =  0       (eps = 0: 厳密ゼロ境界)
     thresh = +epsVal  (eps < 0: より厳しい条件)

   epsVal = 10^-10 で固定し, thresh の符号のみ変化させる.

   全5種の正多面体 x R1/R2/R3 x fallback あり/なし x thresh 3種
   の成功ペア数をまとめて表示する.
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "peeling3DLoxo.m"];

epsVal = 10^-10;

(* ----------------------------------------------------------------
   閾値付き fallback あり版 (thresh パラメータ)
   peeling3DLoxo.m の p3lPeelPair を thresh を引数にして再実装
   ---------------------------------------------------------------- *)
p3lPeelPairThresh[cc_, adj_, nf_, top_, f2_, rule_, thresh_] :=
  Module[{order, inOrder, nbrs, prev, last,
          pending, leftCands, nHat, fHat, lHat,
          angles, absAngles, zvals, maxz, minz, tied, tAngles, next},

    order   = {top, f2};
    inOrder = <|top -> True, f2 -> True|>;
    nbrs    = Table[DeleteCases[DeleteCases[adj[[i]], top], f2], {i, nf}];
    prev = top; last = f2;

    While[Length[order] < nf,
      pending = Select[nbrs[[last]], !KeyExistsQ[inOrder, #]&];
      nbrs[[last]] = pending;
      If[Length[pending] == 0, Break[]];

      {nHat, fHat, lHat} = p3lDarbouxFrame[cc[[prev]], cc[[last]]];

      leftCands = Select[pending, Function[j,
        Det[{cc[[top]], cc[[last]], cc[[j]]}] >= thresh
      ]];

      next = Which[
        rule === "maxphi" && Length[leftCands] > 0,
          angles = p3lAngleFromForward[fHat, lHat, cc[[#]]] & /@ leftCands;
          leftCands[[ First[Ordering[-angles, 1]] ]],

        rule === "loxo" && Length[leftCands] > 0,
          angles = p3lAngleFromForward[fHat, lHat, cc[[#]]] & /@ leftCands;
          leftCands[[ First[Ordering[angles, 1]] ]],

        rule === "maxz" && Length[leftCands] > 0,
          zvals  = cc[[leftCands, 3]];
          maxz   = Max[zvals];
          tied   = Select[leftCands, cc[[#, 3]] >= maxz - 10^-10 &];
          tAngles = p3lAngleFromForward[fHat, lHat, cc[[#]]] & /@ tied;
          tied[[ First[Ordering[-tAngles, 1]] ]],

        (* fallback: 左候補なし -> pending 全体から選択 *)
        rule === "maxphi" && Length[pending] > 0,
          absAngles = Abs[p3lAngleFromForward[fHat, lHat, cc[[#]]]] & /@ pending;
          pending[[ First[Ordering[-absAngles, 1]] ]],

        rule === "loxo" && Length[pending] > 0,
          absAngles = Abs[p3lAngleFromForward[fHat, lHat, cc[[#]]]] & /@ pending;
          pending[[ First[Ordering[absAngles, 1]] ]],

        rule === "maxz" && Length[pending] > 0,
          zvals  = cc[[pending, 3]];
          minz   = Min[zvals];
          tied   = Select[pending, cc[[#, 3]] <= minz + 10^-10 &];
          tAngles = Abs[p3lAngleFromForward[fHat, lHat, cc[[#]]]] & /@ tied;
          tied[[ First[Ordering[-tAngles, 1]] ]],

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

(* ----------------------------------------------------------------
   閾値付き fallback なし版
   ---------------------------------------------------------------- *)
p3lPeelPairNoFBThresh[cc_, adj_, nf_, top_, f2_, rule_, thresh_] :=
  Module[{order, inOrder, nbrs, prev, last,
          pending, leftCands, nHat, fHat, lHat,
          angles, zvals, maxz, tied, tAngles, next},

    order   = {top, f2};
    inOrder = <|top -> True, f2 -> True|>;
    nbrs    = Table[DeleteCases[DeleteCases[adj[[i]], top], f2], {i, nf}];
    prev = top; last = f2;

    While[Length[order] < nf,
      pending = Select[nbrs[[last]], !KeyExistsQ[inOrder, #]&];
      nbrs[[last]] = pending;
      If[Length[pending] == 0, Break[]];

      {nHat, fHat, lHat} = p3lDarbouxFrame[cc[[prev]], cc[[last]]];

      leftCands = Select[pending, Det[{cc[[top]], cc[[last]], cc[[#]]}] >= thresh &];

      If[Length[leftCands] == 0, Break[]];

      next = Which[
        rule === "maxphi",
          angles = p3lAngleFromForward[fHat, lHat, cc[[#]]] & /@ leftCands;
          leftCands[[ First[Ordering[-angles, 1]] ]],
        rule === "loxo",
          angles = p3lAngleFromForward[fHat, lHat, cc[[#]]] & /@ leftCands;
          leftCands[[ First[Ordering[angles, 1]] ]],
        rule === "maxz",
          zvals  = cc[[leftCands, 3]];
          maxz   = Max[zvals];
          tied   = Select[leftCands, cc[[#, 3]] >= maxz - 10^-10 &];
          tAngles = p3lAngleFromForward[fHat, lHat, cc[[#]]] & /@ tied;
          tied[[ First[Ordering[-tAngles, 1]] ]],
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

(* ----------------------------------------------------------------
   全ペア実行 (閾値付き, with fallback)
   ---------------------------------------------------------------- *)
p3lRunAllThresh[vers_, faces_, rule_, thresh_] :=
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
          res = p3lPeelPairThresh[cc, adj, nf, top, f2Cands[[fi]], rule, thresh];
          AppendTo[results, {top, f2Cands[[fi]], res[[1]], res[[2]]}]
        ],
        {fi, Length[f2Cands]}
      ],
      {top, nf}
    ];
    results
  ];

(* ----------------------------------------------------------------
   全ペア実行 (閾値付き, no fallback)
   ---------------------------------------------------------------- *)
p3lRunAllNoFBThresh[vers_, faces_, rule_, thresh_] :=
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
          res = p3lPeelPairNoFBThresh[cc, adj, nf, top, f2Cands[[fi]], rule, thresh];
          AppendTo[results, {top, f2Cands[[fi]], res[[1]], res[[2]]}]
        ],
        {fi, Length[f2Cands]}
      ],
      {top, nf}
    ];
    results
  ];

(* ================================================================
   メインループ
   ================================================================ *)
solidNames = {"Tetrahedron", "Cube", "Octahedron", "Dodecahedron", "Icosahedron"};
rules      = {"maxphi", "loxo", "maxz"};
ruleKeys   = {"R1", "R2", "R3"};

(* 閾値の3ケース *)
threshCases = {
  {-epsVal,  "eps>0 (thresh=-e)", "e+"},
  {0,        "eps=0 (thresh=0) ", "e0"},
  {+epsVal,  "eps<0 (thresh=+e)", "e-"}
};

Print[""];
Print[StringRepeat["=", 90]];
Print["  Platonic Solids: 左条件の閾値 (eps) の影響比較"];
Print["  条件: det(c_1, c_k, c_j) >= thresh"];
Print["  thresh の3ケース:"];
Print["    eps>0: thresh = -10^-10  (現在の設定, 許容あり)"];
Print["    eps=0: thresh =  0       (厳密ゼロ境界)"];
Print["    eps<0: thresh = +10^-10  (より厳しい条件)"];
Print[StringRepeat["=", 90]];

Do[
  rule  = rules[[ri]];
  rkey  = ruleKeys[[ri]];

  Print[""];
  Print[StringRepeat["-", 90]];
  Print["  規則: ", rkey, " (", rule, ")"];
  Print[StringRepeat["-", 90]];
  Print[StringPadLeft["Solid", 14], "  ",
        StringPadLeft["F", 3], "  ",
        StringPadLeft["Pairs", 5], "  ",
        " -- with fallback --        -- no fallback --"];
  Print[StringPadLeft["", 25], "  ",
        StringPadLeft["eps>0", 6], " ",
        StringPadLeft["eps=0", 6], " ",
        StringPadLeft["eps<0", 6], "  ",
        StringPadLeft["eps>0", 6], " ",
        StringPadLeft["eps=0", 6], " ",
        StringPadLeft["eps<0", 6]];
  Print[StringRepeat["-", 90]];

  Do[
    {vers, faces} = p3lExtractPolyhedron[name];
    nf  = Length[faces];
    adj = p3lBuildAdj[faces];
    totalPairs = Total[Length /@ adj];

    wCounts = {};
    nCounts = {};
    Do[
      thresh = threshCases[[ti, 1]];
      rW = p3lRunAllThresh[vers, faces, rule, thresh];
      nW = Length[Select[rW, #[[4]]&]];
      AppendTo[wCounts, nW];

      rN = p3lRunAllNoFBThresh[vers, faces, rule, thresh];
      nN = Length[Select[rN, #[[4]]&]];
      AppendTo[nCounts, nN],
      {ti, 3}
    ];

    Print[StringPadLeft[name, 14], "  ",
          StringPadLeft[ToString[nf], 3], "  ",
          StringPadLeft[ToString[totalPairs], 5], "  ",
          StringPadLeft[ToString[wCounts[[1]]], 6], " ",
          StringPadLeft[ToString[wCounts[[2]]], 6], " ",
          StringPadLeft[ToString[wCounts[[3]]], 6], "  ",
          StringPadLeft[ToString[nCounts[[1]]], 6], " ",
          StringPadLeft[ToString[nCounts[[2]]], 6], " ",
          StringPadLeft[ToString[nCounts[[3]]], 6]],
    {name, solidNames}
  ],
  {ri, 3}
];

Print[""];
Print[StringRepeat["=", 90]];
Print["凡例: 成功ペア数 (/ Pairs で割ると成功率)"];
Print["  eps>0: det >= -10^-10  (許容あり, 現在の設定)"];
Print["  eps=0: det >=  0       (厳密ゼロ境界)"];
Print["  eps<0: det >= +10^-10  (より厳しい)"];
Print[""];
Print["Done."];
