(* ================================================================
   run_platonic_fallback.m
   正5種多面体 × R1/R2/R3 について
   fallback あり版（p3lRunAll）と fallback なし版を比較する．
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "peeling3DLoxo.m"];

(* --- fallback なし版 1ペア関数 --- *)
p3lPeelPairNoFB[cc_, adj_, nf_, top_, f2_, rule_] :=
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
      (* 候補が1つだけの場合は左条件によらず選択する *)
      leftCands = If[Length[pending] == 1,
        pending,
        Select[pending, Det[{cc[[prev]], cc[[last]], cc[[#]]}] >= 0 &]
      ];

      (* fallback なし: 左候補なければ即終了 *)
      If[Length[leftCands] == 0, Break[]];

      next = Which[
        rule === "maxphi",
          angles = p3lAngleFromForward[fHat, lHat, cc[[#]]] & /@ leftCands;
          leftCands[[ First[Ordering[-angles, 1]] ]],
        rule === "loxo",
          angles = p3lAngleFromForward[fHat, lHat, cc[[#]]] & /@ leftCands;
          leftCands[[ First[Ordering[angles,  1]] ]],
        rule === "maxz",
          zvals   = cc[[leftCands, 3]];
          maxz    = Max[zvals];
          tied    = Select[leftCands, cc[[#, 3]] >= maxz - 10^-10 &];
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

(* --- fallback なし版 全ペア実行 --- *)
p3lRunAllNoFB[vers_, faces_, rule_] :=
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
          res = p3lPeelPairNoFB[cc, adj, nf, top, f2Cands[[fi]], rule];
          AppendTo[results, {top, f2Cands[[fi]], res[[1]], res[[2]]}]
        ],
        {fi, Length[f2Cands]}
      ],
      {top, nf}
    ];
    results
  ];

(* --- 正多面体5種 --- *)
solidNames = {"Tetrahedron", "Cube", "Octahedron", "Dodecahedron", "Icosahedron"};
rules      = {"maxphi", "loxo", "maxz"};
ruleLabels = {"R1", "R2", "R3"};

Print[""];
Print[StringRepeat["=", 72]];
Print["  Platonic Solids: fallback あり vs なし"];
Print[StringRepeat["=", 72]];
Print[StringPadLeft["Solid", 13], "  ",
      StringPadLeft["F",  3], "  ",
      StringPadLeft["Pairs", 5], "  ",
      "  R1(w)  R1(n)  R2(w)  R2(n)  R3(w)  R3(n)"];
Print[StringRepeat["-", 72]];

allRes = <||>;

Do[
  {vers, faces} = p3lExtractPolyhedron[name];
  nf  = Length[faces];
  adj = p3lBuildAdj[faces];
  totalPairs = Total[Length /@ adj];

  row = {name, nf, totalPairs};

  Do[
    (* with fallback *)
    rWith = p3lRunAll[vers, faces, rules[[ri]]];
    nWith = Length[Select[rWith, #[[4]]&]];
    (* no fallback *)
    rNoFB = p3lRunAllNoFB[vers, faces, rules[[ri]]];
    nNoFB = Length[Select[rNoFB, #[[4]]&]];
    AppendTo[row, {nWith, nNoFB}],
    {ri, 3}
  ];

  allRes[name] = row;

  Print[StringPadLeft[name, 13], "  ",
        StringPadLeft[ToString[nf], 3], "  ",
        StringPadLeft[ToString[totalPairs], 5], "  ",
        Sequence @@ Table[
          StringPadLeft[ToString[row[[3+ri,1]]], 6] <>
          StringPadLeft[ToString[row[[3+ri,2]]], 6],
          {ri, 3}
        ]
  ],
  {name, solidNames}
];

Print[StringRepeat["-", 72]];
Print["(w) = with fallback,  (n) = no fallback"];
Print["\nDone."];
