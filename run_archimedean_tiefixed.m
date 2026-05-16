(* ================================================================
   run_archimedean_tiefixed.m
   ================================================================
   インデックス依存タイブレーク修正後のコードでアルキメデス13種を
   with-fallback / no-fallback の両条件で再計算し，
   旧結果と比較する．

   変更点（peeling3DLoxo.m および run_archimedean_nofallback.m）:
     RS: max φ → φ タイは max Det
     R2: min φ → φ タイは min Det
     RZ: max z → max φ → max Det の三段階
     各フォールバック路も同様に Det 基準を追加

   出力:
     archimedean_faceup_v2.mx      -- with-fallback 新結果
     archimedean_nofallback_v2.mx  -- no-fallback  新結果
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];

(* peeling3DLoxo.m（修正済み）だけをロード．
   run_archimedean_nofallback.m は Get しない（旧 .mx を上書きするため）． *)
Get[scriptDir <> "peeling3DLoxo.m"];

(* ----------------------------------------------------------------
   no-fallback 版の関数定義（run_archimedean_nofallback.m の修正済み版と同一）
   ---------------------------------------------------------------- *)
p3lPeelPairNoFB2[cc_, adj_, nf_, top_, f2_, rule_] :=
  Module[{order, inOrder, nbrs, prev, last,
          pending, leftCands, nHat, fHat, lHat,
          angles, zvals, maxz, maxPhi,
          tied, tiedPhi, tAngles, tDets, next,
          eps},
    eps = 10^-10;

    order   = {top, f2};
    inOrder = <|top -> True, f2 -> True|>;
    nbrs    = Table[
      DeleteCases[DeleteCases[adj[[i]], top], f2],
      {i, nf}
    ];
    prev = top;
    last = f2;

    While[Length[order] < nf,
      pending = Select[nbrs[[last]], !KeyExistsQ[inOrder, #] &];
      nbrs[[last]] = pending;
      If[Length[pending] == 0, Break[]];

      {nHat, fHat, lHat} = p3lDarbouxFrame[cc[[prev]], cc[[last]]];

      leftCands = Select[pending, Function[j,
        Det[{cc[[top]], cc[[last]], cc[[j]]}] >= -eps
      ]];

      If[Length[leftCands] == 0, Break[]];

      next = Which[
        rule === "maxphi",
          angles = p3lAngleFromForward[fHat, lHat, cc[[#]]] & /@ leftCands;
          maxPhi = Max[angles];
          tied   = Pick[leftCands, Thread[angles >= maxPhi - 10^-10]];
          If[Length[tied] == 1, tied[[1]],
            tDets = Det[{cc[[top]], cc[[last]], cc[[#]]}] & /@ tied;
            tied[[ First[Ordering[-tDets, 1]] ]]
          ],

        rule === "maxz",
          zvals = cc[[leftCands, 3]];
          maxz  = Max[zvals];
          tied  = Select[leftCands, cc[[#, 3]] >= maxz - 10^-10 &];
          If[Length[tied] == 1, tied[[1]],
            tAngles = p3lAngleFromForward[fHat, lHat, cc[[#]]] & /@ tied;
            maxPhi  = Max[tAngles];
            tiedPhi = Pick[tied, Thread[tAngles >= maxPhi - 10^-10]];
            If[Length[tiedPhi] == 1, tiedPhi[[1]],
              tDets = Det[{cc[[top]], cc[[last]], cc[[#]]}] & /@ tiedPhi;
              tiedPhi[[ First[Ordering[-tDets, 1]] ]]
            ]
          ],

        True, Break[]; Null
      ];

      If[next === Null, Break[]];
      AppendTo[order, next];
      inOrder[next] = True;
      nbrs = Table[DeleteCases[nbrs[[i]], next], {i, nf}];
      prev = last;
      last = next;
    ];
    {order, Length[order] == nf}
  ];

p3lRunAllNoFB2[vers_, faces_, rule_] :=
  Module[{nf, adj, lv, cc, results, f2Cands, order, complete},
    nf      = Length[faces];
    adj     = p3lBuildAdj[faces];
    results = {};
    Do[
      lv      = p3lAlignTopToZ[vers, faces, top];
      cc      = p3lFaceCentroid[#, lv, faces] & /@ Range[nf];
      f2Cands = adj[[top]];
      Do[
        {order, complete} =
          p3lPeelPairNoFB2[cc, adj, nf, top, f2Cands[[f2i]], rule];
        AppendTo[results, {top, f2i, order, complete}],
        {f2i, Length[f2Cands]}
      ],
      {top, nf}
    ];
    results
  ];

(* ----------------------------------------------------------------
   旧結果の読み込み
   ---------------------------------------------------------------- *)
oldWF = <||>;
oldNF = <||>;

wfFile = scriptDir <> "archimedean_faceup_results.mx";
nfFile = scriptDir <> "archimedean_nofallback_results.mx";

If[FileExistsQ[wfFile],
  Get[wfFile]; oldWF = allResults;
  Print["Loaded old with-FB: ", wfFile],
  Print["WARNING: archimedean_faceup_results.mx not found"]
];
If[FileExistsQ[nfFile],
  Get[nfFile]; oldNF = nfResults;
  Print["Loaded old no-FB:   ", nfFile],
  Print["WARNING: archimedean_nofallback_results.mx not found"]
];

(* ----------------------------------------------------------------
   アルキメデス13種の再計算（新コード）
   ---------------------------------------------------------------- *)
solidNames = {
  "TruncatedTetrahedron", "Cuboctahedron", "TruncatedCube",
  "TruncatedOctahedron", "Rhombicuboctahedron", "TruncatedCuboctahedron",
  "SnubCube", "Icosidodecahedron", "TruncatedDodecahedron",
  "TruncatedIcosahedron", "Rhombicosidodecahedron",
  "TruncatedIcosidodecahedron", "SnubDodecahedron"
};

Print[""];
Print[StringRepeat["=", 70]];
Print["  Re-computing with tie-fixed code ..."];
Print[StringRepeat["=", 70]];

newWF = <||>;
newNF = <||>;

Do[
  {vers, faces} = p3lExtractPolyhedron[name];
  nf  = Length[faces];
  adj = p3lBuildAdj[faces];
  totalPairs = Total[Length /@ adj];

  Print["\n-- ", name, "  (F=", nf, ", pairs=", totalPairs, ") --"];

  (* with-fallback: p3lRunAll（peeling3DLoxo.m の修正済み版） *)
  r1w = p3lRunAll[vers, faces, "maxphi"];
  r2w = p3lRunAll[vers, faces, "loxo"];
  r3w = p3lRunAll[vers, faces, "maxz"];
  n1w = Length[Select[r1w, #[[4]] &]];
  n2w = Length[Select[r2w, #[[4]] &]];
  n3w = Length[Select[r3w, #[[4]] &]];
  Print["  with-FB  RS:", n1w, "  R2:", n2w, "  RZ:", n3w];

  (* no-fallback: p3lRunAllNoFB2（本スクリプトで定義した修正済み版） *)
  r1n = p3lRunAllNoFB2[vers, faces, "maxphi"];
  r3n = p3lRunAllNoFB2[vers, faces, "maxz"];
  n1n = Length[Select[r1n, #[[4]] &]];
  n3n = Length[Select[r3n, #[[4]] &]];
  Print["  no-FB    RS:", n1n, "          RZ:", n3n];

  newWF[name] = <|
    "nFaces" -> nf, "totalPairs" -> totalPairs,
    "R1" -> n1w, "R2" -> n2w, "R3" -> n3w,
    "details" -> <|"R1" -> r1w, "R2" -> r2w, "R3" -> r3w|>
  |>;
  newNF[name] = <|
    "nFaces" -> nf, "totalPairs" -> totalPairs,
    "RS_noFB" -> n1n, "RZ_noFB" -> n3n,
    "details" -> <|"RS" -> r1n, "RZ" -> r3n|>
  |>;,
  {name, solidNames}
];

(* ----------------------------------------------------------------
   結果の保存（旧 .mx は上書きしない）
   ---------------------------------------------------------------- *)
allResults = newWF;
DumpSave[scriptDir <> "archimedean_faceup_v2.mx", allResults];
Print["\nSaved: archimedean_faceup_v2.mx"];

nfResults = newNF;
DumpSave[scriptDir <> "archimedean_nofallback_v2.mx", nfResults];
Print["Saved: archimedean_nofallback_v2.mx"];

(* ----------------------------------------------------------------
   旧 vs 新 比較表
   ---------------------------------------------------------------- *)
Print[""];
Print[StringRepeat["=", 105]];
Print["  Comparison: old vs tie-fixed (new)   (*) = changed"];
Print[StringRepeat["=", 105]];
Print[StringPadLeft["Solid", 30], "  Pairs",
      "  RS(w)old->new  R2(w)old->new  RZ(w)old->new",
      "  RS(n)old->new  RZ(n)old->new"];
Print[StringRepeat["-", 105]];

anyChange = False;

fmtCell[old_, new_] :=
  StringPadLeft[
    StringJoin[
      If[old === "?", "?", ToString[old]], "->", ToString[new],
      If[old =!= "?" && old =!= new, "(*)", "   "]
    ], 13
  ];

Do[
  tp   = newWF[name]["totalPairs"];
  n1wN = newWF[name]["R1"];
  n2wN = newWF[name]["R2"];
  n3wN = newWF[name]["R3"];
  n1nN = newNF[name]["RS_noFB"];
  n3nN = newNF[name]["RZ_noFB"];

  {n1wO, n2wO, n3wO} = If[KeyExistsQ[oldWF, name],
    {oldWF[name]["R1"], oldWF[name]["R2"], oldWF[name]["R3"]},
    {"?", "?", "?"}
  ];
  {n1nO, n3nO} = If[KeyExistsQ[oldNF, name],
    {oldNF[name]["RS_noFB"], oldNF[name]["RZ_noFB"]},
    {"?", "?"}
  ];

  changed = AnyTrue[
    {n1wO =!= "?" && n1wO =!= n1wN,
     n2wO =!= "?" && n2wO =!= n2wN,
     n3wO =!= "?" && n3wO =!= n3wN,
     n1nO =!= "?" && n1nO =!= n1nN,
     n3nO =!= "?" && n3nO =!= n3nN},
    TrueQ
  ];
  If[changed, anyChange = True];

  Print[
    StringPadLeft[name, 30], "  ", StringPadLeft[ToString[tp], 5], "  ",
    fmtCell[n1wO, n1wN], "  ", fmtCell[n2wO, n2wN], "  ", fmtCell[n3wO, n3wN], "  ",
    fmtCell[n1nO, n1nN], "  ", fmtCell[n3nO, n3nN]
  ],
  {name, solidNames}
];

Print[StringRepeat["-", 105]];
If[anyChange,
  Print["  (*) = result changed by tie-fixed code"],
  Print["  No changes -- tie-fixed code gives identical results for all solids."]
];

Print["\nDone."];
