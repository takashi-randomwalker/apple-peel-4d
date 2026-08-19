(* ================================================================
   run_platonic_updated.m
   R3 z タイブレーク修正後の正多面体 R1/R2/R3 × fallback あり/なし の再計算
   結果は dataPlatonic.mx に保存する．

   保存データ構造:
     platonicData[name] = <|
       "nFaces"      -> 面数,
       "totalPairs"  -> (top,f2) ペアの総数,
       "R1" -> <| "withFallback" -> {{top, f2, order, success}, ...},
                  "noFallback"   -> {{top, f2, order, success}, ...} |>,
       "R2" -> <| ... |>,
       "R3" -> <| ... |>
     |>

   各エントリの形式:
     top     : 開始面のインデックス (Integer)
     f2      : 2番目の面のインデックス (Integer, 面番号そのもの)
     order   : ペーリング順序 (List of Integer, 成功時は全面数分, 失敗時は途中まで)
     success : 完全ペーリングに成功したか (True / False)
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "peeling3DLoxo.m"];

(* ----------------------------------------------------------------
   fallback なし版（単一候補優先・R3 max φ タイブレークを適用）
   戻り値の形式: {top, f2実面番号, order, success}
   ---------------------------------------------------------------- *)
p3lPeelPairNoFB2[cc_, adj_, nf_, top_, f2_, rule_] :=
  Module[{order, inOrder, nbrs, last,
          pending, leftCands, dets, zvals, maxz, tied, tDets, next, eps},
    eps = 10^-10;

    order   = {top, f2};
    inOrder = <|top -> True, f2 -> True|>;
    nbrs    = Table[DeleteCases[DeleteCases[adj[[i]], top], f2], {i, nf}];
    last = f2;

    While[Length[order] < nf,
      pending = Select[nbrs[[last]], !KeyExistsQ[inOrder, #]&];
      nbrs[[last]] = pending;
      If[Length[pending] == 0, Break[]];

      leftCands = Select[pending, Det[{cc[[top]], cc[[last]], cc[[#]]}] <= eps &];

      If[Length[leftCands] == 0, Break[]];

      next = Which[
        rule === "maxphi",
          dets = Det[{cc[[top]], cc[[last]], cc[[#]]}] & /@ leftCands;
          leftCands[[ First[Ordering[dets, 1]] ]],
        rule === "loxo",
          dets = Det[{cc[[top]], cc[[last]], cc[[#]]}] & /@ leftCands;
          leftCands[[ First[Ordering[-dets, 1]] ]],
        rule === "maxz",
          zvals = cc[[leftCands, 3]];
          maxz  = Max[zvals];
          tied  = Select[leftCands, cc[[#, 3]] >= maxz - 10^-10 &];
          If[Length[tied] == 1, tied[[1]],
            tDets = Det[{cc[[top]], cc[[last]], cc[[#]]}] & /@ tied;
            tied[[ First[Ordering[tDets, 1]] ]]
          ],
        True, Break[]; Null
      ];
      If[next === Null, Break[]];
      AppendTo[order, next];
      inOrder[next] = True;
      nbrs = Table[DeleteCases[nbrs[[i]], next], {i, nf}];
      last = next;
    ];
    {order, Length[order] == nf}
  ];

p3lRunAllNoFB2[vers_, faces_, rule_] :=
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
          res = p3lPeelPairNoFB2[cc, adj, nf, top, f2Cands[[fi]], rule];
          (* f2 は実際の面番号 f2Cands[[fi]] を格納 *)
          AppendTo[results, {top, f2Cands[[fi]], res[[1]], res[[2]]}]
        ],
        {fi, Length[f2Cands]}
      ],
      {top, nf}
    ];
    results
  ];

(* ----------------------------------------------------------------
   p3lRunAll の戻り値は {top, f2i, order, complete} で
   f2i は adj[[top]] 内のインデックスなので，実面番号に変換する
   ---------------------------------------------------------------- *)
p3lRunAllFixed[vers_, faces_, rule_] :=
  Module[{adj, raw},
    adj = p3lBuildAdj[faces];
    raw = p3lRunAll[vers, faces, rule];
    (* {top, f2i, order, complete} -> {top, f2, order, complete} *)
    {#[[1]], adj[[ #[[1]] ]][[ #[[2]] ]], #[[3]], #[[4]]} & /@ raw
  ];

(* ----------------------------------------------------------------
   正多面体 5 種
   ---------------------------------------------------------------- *)
solidNames = {"Tetrahedron", "Cube", "Octahedron", "Dodecahedron", "Icosahedron"};
rules      = {"maxphi", "loxo", "maxz"};
ruleKeys   = {"R1", "R2", "R3"};

Print[""];
Print[StringRepeat["=", 72]];
Print["  Platonic Solids: R3 タイブレーク修正後の結果"];
Print[StringRepeat["=", 72]];
Print[StringPadLeft["Solid", 13], "  ",
      StringPadLeft["F",  3], "  ",
      StringPadLeft["Pairs", 5], "  ",
      "  R1(w)  R1(n)  R2(w)  R2(n)  R3(w)  R3(n)"];
Print[StringRepeat["-", 72]];

platonicData = <||>;

Do[
  {vers, faces} = p3lExtractPolyhedron[name];
  nf  = Length[faces];
  adj = p3lBuildAdj[faces];
  totalPairs = Total[Length /@ adj];

  row     = {name, nf, totalPairs};
  solidEntry = <|"nFaces" -> nf, "totalPairs" -> totalPairs|>;

  Do[
    rWith = p3lRunAllFixed[vers, faces, rules[[ri]]];
    nWith = Length[Select[rWith, #[[4]]&]];

    rNoFB = p3lRunAllNoFB2[vers, faces, rules[[ri]]];
    nNoFB = Length[Select[rNoFB, #[[4]]&]];

    AppendTo[row, {nWith, nNoFB}];

    solidEntry[ruleKeys[[ri]]] = <|
      "withFallback" -> rWith,
      "noFallback"   -> rNoFB
    |>,
    {ri, 3}
  ];

  platonicData[name] = solidEntry;

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
Print["単一候補優先・R3 z タイブレーク(max φ = R1)は全ケースに適用"];

(* ----------------------------------------------------------------
   結果の保存
   ---------------------------------------------------------------- *)
outFile = scriptDir <> "dataPlatonic.mx";
DumpSave[outFile, platonicData];
Print["\nSaved: ", outFile];
Print["変数名: platonicData"];
Print["  platonicData[name][rule][\"withFallback\" | \"noFallback\"]"];
Print["  各エントリ: {top, f2, order, success}"];
Print["  order は成功時は全面，失敗時は途中まで"];
Print["\nDone."];
