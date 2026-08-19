(* ================================================================
   run_archimedean_nofallback.m
   ================================================================
   アルキメデス13種を **フォールバックなし** 条件で計算し，
   RS (maxphi) / RZ (maxz) の成功数を集計する．
   既存の with-fallback 結果 (archimedean_faceup_results.mx) と比較する．
   結果は archimedean_nofallback_results.mx に保存する．
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "peeling3DLoxo.m"];

(* ----------------------------------------------------------------
   フォールバックなし版の peelPair
   leftCands が空のとき即座に失敗（Break）する．
   左条件・eps・グローバル c_1 参照はすべて with-fallback 版と同じ．
   ---------------------------------------------------------------- *)
p3lPeelPairNoFB[cc_, adj_, nf_, top_, f2_, rule_] :=
  Module[{order, inOrder, nbrs, last,
          pending, leftCands, dets, zvals, maxz,
          tied, tDets, next,
          eps},
    eps = 10^-10;

    order   = {top, f2};
    inOrder = <|top -> True, f2 -> True|>;
    nbrs    = Table[
      DeleteCases[DeleteCases[adj[[i]], top], f2],
      {i, nf}
    ];
    last = f2;

    While[Length[order] < nf,

      pending = Select[nbrs[[last]], !KeyExistsQ[inOrder, #]&];
      nbrs[[last]] = pending;
      If[Length[pending] == 0, Break[]];

      leftCands = Select[pending, Function[j,
        Det[{cc[[top]], cc[[last]], cc[[j]]}] <= eps
      ]];

      (* フォールバックなし: leftCands が空なら即失敗 *)
      If[Length[leftCands] == 0, Break[]];

      next = Which[
        rule === "maxphi",
          (* RS: min Det（最も負 = 最大右旋回） *)
          dets = Det[{cc[[top]], cc[[last]], cc[[#]]}] & /@ leftCands;
          leftCands[[ First[Ordering[dets, 1]] ]],

        rule === "maxz",
          (* RZ: max z → min Det タイブレーク *)
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

p3lRunAllNoFB[vers_, faces_, rule_] :=
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
          p3lPeelPairNoFB[cc, adj, nf, top, f2Cands[[f2i]], rule];
        AppendTo[results, {top, f2i, order, complete}],
        {f2i, Length[f2Cands]}
      ],
      {top, nf}
    ];
    results
  ];


(* --- アルキメデス13種 --- *)
solidNames = {
  "TruncatedTetrahedron",
  "Cuboctahedron",
  "TruncatedCube",
  "TruncatedOctahedron",
  "Rhombicuboctahedron",
  "TruncatedCuboctahedron",
  "SnubCube",
  "Icosidodecahedron",
  "TruncatedDodecahedron",
  "TruncatedIcosahedron",
  "Rhombicosidodecahedron",
  "TruncatedIcosidodecahedron",
  "SnubDodecahedron"
};

(* 既存の with-fallback 結果を読み込む *)
wfFile = scriptDir <> "archimedean_faceup_results.mx";
If[FileExistsQ[wfFile],
  Get[wfFile];
  hasWF = True,
  hasWF = False;
  Print["WARNING: archimedean_faceup_results.mx not found; skipping comparison."]
];

Print[""];
Print[StringRepeat["=", 70]];
Print["  Archimedean Solids: No-Fallback (RS and RZ)"];
Print[StringRepeat["=", 70]];

nfResults = <||>;

Do[
  {vers, faces} = p3lExtractPolyhedron[name];
  nf  = Length[faces];
  adj = p3lBuildAdj[faces];
  totalPairs = Total[Length /@ adj];

  Print["\n-- ", name, "  (F=", nf, ", pairs=", totalPairs, ") --"];

  solidRes = <||>;

  (* RS (maxphi) *)
  t0 = AbsoluteTime[];
  rs = p3lRunAllNoFB[vers, faces, "maxphi"];
  dt = AbsoluteTime[] - t0;
  nRS = Length[Select[rs, #[[4]]&]];
  Print["  RS (no-FB): ", nRS, "/", totalPairs,
        "  (", NumberForm[N[100 nRS/totalPairs], {4,1}], "%)",
        "  ", NumberForm[dt, {4,2}], " s"];
  solidRes["RS"] = rs;

  (* RZ (maxz) *)
  t0 = AbsoluteTime[];
  rz = p3lRunAllNoFB[vers, faces, "maxz"];
  dt = AbsoluteTime[] - t0;
  nRZ = Length[Select[rz, #[[4]]&]];
  Print["  RZ (no-FB): ", nRZ, "/", totalPairs,
        "  (", NumberForm[N[100 nRZ/totalPairs], {4,1}], "%)",
        "  ", NumberForm[dt, {4,2}], " s"];
  solidRes["RZ"] = rz;

  nfResults[name] = <|
    "nFaces" -> nf, "totalPairs" -> totalPairs,
    "RS_noFB" -> nRS, "RZ_noFB" -> nRZ,
    "details" -> solidRes
  |>;,
  {name, solidNames}
];

(* --- 結果の保存 --- *)
outFile = scriptDir <> "archimedean_nofallback_results.mx";
DumpSave[outFile, nfResults];
Print["\nSaved: ", outFile];

(* --- サマリ比較表 --- *)
Print[""];
Print[StringRepeat["=", 80]];
Print["  Summary: With-Fallback vs No-Fallback"];
Print[StringRepeat["=", 80]];
Print[StringPadLeft["Solid", 30], "  ",
      StringPadLeft["Pairs", 5], "  ",
      StringPadLeft["RS(w)", 6], "  ",
      StringPadLeft["RS(n)", 6], "  ",
      StringPadLeft["RZ(w)", 6], "  ",
      StringPadLeft["RZ(n)", 6]];
Print[StringRepeat["-", 65]];
Do[
  r  = nfResults[name];
  tp = r["totalPairs"];
  nRS = r["RS_noFB"];
  nRZ = r["RZ_noFB"];
  If[hasWF,
    wRS = allResults[name]["R1"];  (* with-FB の R1 = RS *)
    wRZ = allResults[name]["R3"];  (* with-FB の R3 = RZ *)
    Print[StringPadLeft[name, 30], "  ",
          StringPadLeft[ToString[tp], 5], "  ",
          StringPadLeft[ToString[wRS], 6], "  ",
          StringPadLeft[ToString[nRS], 6], "  ",
          StringPadLeft[ToString[wRZ], 6], "  ",
          StringPadLeft[ToString[nRZ], 6]],
    Print[StringPadLeft[name, 30], "  ",
          StringPadLeft[ToString[tp], 5], "  ",
          StringPadLeft["n/a", 6], "  ",
          StringPadLeft[ToString[nRS], 6], "  ",
          StringPadLeft["n/a", 6], "  ",
          StringPadLeft[ToString[nRZ], 6]]
  ],
  {name, solidNames}
];
Print[StringRepeat["-", 65]];
Print["  w = with fallback, n = no fallback"];

Print["\nDone."];
