(* ================================================================
   run_platonic_strict.m
   厳密右条件（eps < 0, thresh = -10^-10）での正多面体 RS/RZ 計算
   右条件: det(c_1, c_k, c_j) <= -epsVal

   保存先: dataPlatonic_strict.mx
   変数名: platonicDataStrict

   データ構造（dataPlatonic.mx と同形式）:
     platonicDataStrict[name] = <|
       "nFaces"      -> 面数,
       "totalPairs"  -> ペア総数,
       "RS" -> <| "withFallback" -> {{top,f2,order,success},...},
                  "noFallback"   -> {{top,f2,order,success},...} |>,
       "RZ" -> <| ... |>
     |>
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "peeling3DLoxo.m"];

epsVal = 10^-10;
thresh = -epsVal;   (* 厳密右条件: det <= -10^-10 *)

(* ----------------------------------------------------------------
   厳密条件 with fallback 版
   ---------------------------------------------------------------- *)
p3lPeelPairStrict[cc_, adj_, nf_, top_, f2_, rule_] :=
  Module[{order, inOrder, nbrs, prev, last,
          pending, leftCands, nHat, fHat, lHat,
          dets, angles, zvals, maxz, minz, tied, tDets, next},

    order   = {top, f2};
    inOrder = <|top -> True, f2 -> True|>;
    nbrs    = Table[DeleteCases[DeleteCases[adj[[i]], top], f2], {i, nf}];
    prev = top; last = f2;

    While[Length[order] < nf,
      pending = Select[nbrs[[last]], !KeyExistsQ[inOrder, #]&];
      nbrs[[last]] = pending;
      If[Length[pending] == 0, Break[]];

      (* 厳密右候補: det <= -epsVal *)
      leftCands = Select[pending, Function[j,
        Det[{cc[[top]], cc[[last]], cc[[j]]}] <= thresh
      ]];

      next = Which[
        rule === "maxphi" && Length[leftCands] > 0,
          dets = Det[{cc[[top]], cc[[last]], cc[[#]]}] & /@ leftCands;
          leftCands[[ First[Ordering[dets, 1]] ]],

        rule === "maxz" && Length[leftCands] > 0,
          zvals = cc[[leftCands, 3]];
          maxz  = Max[zvals];
          tied  = Select[leftCands, cc[[#, 3]] >= maxz - 10^-10 &];
          If[Length[tied] == 1, tied[[1]],
            tDets = Det[{cc[[top]], cc[[last]], cc[[#]]}] & /@ tied;
            tied[[ First[Ordering[tDets, 1]] ]]
          ],

        (* fallback: det > -epsVal の全候補から選択 *)
        rule === "maxphi" && Length[pending] > 0,
          {nHat, fHat, lHat} = p3lDarbouxFrame[cc[[prev]], cc[[last]]];
          angles = p3lAngleFromForward[fHat, lHat, cc[[#]]] & /@ pending;
          pending[[ First[Ordering[angles, 1]] ]],

        rule === "maxz" && Length[pending] > 0,
          zvals  = cc[[pending, 3]];
          minz   = Min[zvals];
          tied   = Select[pending, cc[[#, 3]] <= minz + 10^-10 &];
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
      prev = last; last = next;
    ];
    {order, Length[order] == nf}
  ];

(* ----------------------------------------------------------------
   厳密条件 no fallback 版
   ---------------------------------------------------------------- *)
p3lPeelPairStrictNoFB[cc_, adj_, nf_, top_, f2_, rule_] :=
  Module[{order, inOrder, nbrs, last,
          pending, leftCands, dets, zvals, maxz, tied, tDets, next},

    order   = {top, f2};
    inOrder = <|top -> True, f2 -> True|>;
    nbrs    = Table[DeleteCases[DeleteCases[adj[[i]], top], f2], {i, nf}];
    last = f2;

    While[Length[order] < nf,
      pending = Select[nbrs[[last]], !KeyExistsQ[inOrder, #]&];
      nbrs[[last]] = pending;
      If[Length[pending] == 0, Break[]];

      leftCands = Select[pending, Det[{cc[[top]], cc[[last]], cc[[#]]}] <= thresh &];

      If[Length[leftCands] == 0, Break[]];

      next = Which[
        rule === "maxphi",
          dets = Det[{cc[[top]], cc[[last]], cc[[#]]}] & /@ leftCands;
          leftCands[[ First[Ordering[dets, 1]] ]],
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

(* ----------------------------------------------------------------
   全ペア実行
   ---------------------------------------------------------------- *)
p3lRunAllStrict[vers_, faces_, rule_] :=
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
          res = p3lPeelPairStrict[cc, adj, nf, top, f2Cands[[fi]], rule];
          AppendTo[results, {top, f2Cands[[fi]], res[[1]], res[[2]]}]
        ],
        {fi, Length[f2Cands]}
      ],
      {top, nf}
    ];
    results
  ];

p3lRunAllStrictNoFB[vers_, faces_, rule_] :=
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
          res = p3lPeelPairStrictNoFB[cc, adj, nf, top, f2Cands[[fi]], rule];
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
ruleMap    = <|"RS" -> "maxphi", "RZ" -> "maxz"|>;

Print[""];
Print[StringRepeat["=", 72]];
Print["  Platonic Solids: 厳密右条件 (thresh = -10^-10) の結果"];
Print["  RS = Spiral rule (min phi),  RZ = Zonal rule (max z)"];
Print[StringRepeat["=", 72]];
Print[StringPadLeft["Solid", 13], "  ",
      StringPadLeft["F", 3], "  ",
      StringPadLeft["Pairs", 5], "  ",
      " RS(w)  RS(n)  RZ(w)  RZ(n)"];
Print[StringRepeat["-", 72]];

platonicDataStrict = <||>;

Do[
  {vers, faces} = p3lExtractPolyhedron[name];
  nf  = Length[faces];
  adj = p3lBuildAdj[faces];
  totalPairs = Total[Length /@ adj];

  solidEntry = <|"nFaces" -> nf, "totalPairs" -> totalPairs|>;
  counts = {};

  Do[
    mRule = ruleMap[rkey];

    rWith  = p3lRunAllStrict[vers, faces, mRule];
    nWith  = Length[Select[rWith, #[[4]]&]];

    rNoFB  = p3lRunAllStrictNoFB[vers, faces, mRule];
    nNoFB  = Length[Select[rNoFB, #[[4]]&]];

    AppendTo[counts, nWith];
    AppendTo[counts, nNoFB];

    solidEntry[rkey] = <|
      "withFallback" -> rWith,
      "noFallback"   -> rNoFB
    |>,
    {rkey, {"RS", "RZ"}}
  ];

  platonicDataStrict[name] = solidEntry;

  Print[StringPadLeft[name, 13], "  ",
        StringPadLeft[ToString[nf], 3], "  ",
        StringPadLeft[ToString[totalPairs], 5], "  ",
        StringPadLeft[ToString[counts[[1]]], 6],
        StringPadLeft[ToString[counts[[2]]], 6],
        StringPadLeft[ToString[counts[[3]]], 6],
        StringPadLeft[ToString[counts[[4]]], 6]],
  {name, solidNames}
];

Print[StringRepeat["-", 72]];
Print["(w) = with fallback,  (n) = no fallback"];
Print["thresh = -10^-10: det > -10^-10 の面（境界面含む）は右候補から除外"];

(* 参考: eps > 0 (現行) との比較 *)
Print[""];
Print["参考: 現行 (thresh = -10^-10, eps > 0) との比較"];
Print[StringPadLeft["Solid", 13], "  ",
      StringPadLeft["F", 3], "  ",
      StringPadLeft["Pairs", 5], "  ",
      " RS(w)  RS(n)  RZ(w)  RZ(n)"];
Print[StringRepeat["-", 72]];
Get[scriptDir <> "dataPlatonic.mx"];
Do[
  nf = platonicData[name]["nFaces"];
  tp = platonicData[name]["totalPairs"];
  nRSw = Length[Select[platonicData[name]["R1"]["withFallback"], #[[4]]&]];
  nRSn = Length[Select[platonicData[name]["R1"]["noFallback"],   #[[4]]&]];
  nRZw = Length[Select[platonicData[name]["R3"]["withFallback"], #[[4]]&]];
  nRZn = Length[Select[platonicData[name]["R3"]["noFallback"],   #[[4]]&]];
  Print[StringPadLeft[name, 13], "  ",
        StringPadLeft[ToString[nf], 3], "  ",
        StringPadLeft[ToString[tp], 5], "  ",
        StringPadLeft[ToString[nRSw], 6],
        StringPadLeft[ToString[nRSn], 6],
        StringPadLeft[ToString[nRZw], 6],
        StringPadLeft[ToString[nRZn], 6]],
  {name, solidNames}
];

(* 保存 *)
outFile = scriptDir <> "dataPlatonic_strict.mx";
DumpSave[outFile, platonicDataStrict];
Print[""];
Print["Saved: ", outFile];
Print["変数名: platonicDataStrict"];
Print["  platonicDataStrict[name][\"RS\" | \"RZ\"][\"withFallback\" | \"noFallback\"]"];
Print["  各エントリ: {top, f2, order, success}"];
Print["\nDone."];
