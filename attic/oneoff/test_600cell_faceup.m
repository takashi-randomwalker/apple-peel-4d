(* ================================================================
   test_600cell_faceup.m
   600胞体に 3D-face-centroid-up（共有面重心を +w）+ RZ を適用し，
   全 2400 (C1,C2) ペアの成功/失敗と終了ステップを集計する．
   比較：現行 cell-centroid-up では全 2400 ペアが Impossible．
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];

dataDir = scriptDir <> "_4DData/";
raw   = Get[dataDir <> "f600.m"];
vers  = N[raw[[1]]]; faces = raw[[4]]; cells = raw[[6]];
nc    = Length[cells];
Print["600-cell: ", nc, " cells, ", Length[vers], " vertices, ",
      Length[faces], " faces"];

wHat = {0., 0., 0., 1.};
eps  = 10^-10;

cellVerts[ci_] := Union @@ faces[[ cells[[ci]] ]];

(* 全セル重心 *)
cc0 = Table[Mean[vers[[ cellVerts[ci] ]]], {ci, nc}];

(* face-up 回転行列（方向ベクトル → +w） *)
faceUpRot[dir_] := Module[{d = Normalize[dir - Mean[vers]], ang},
  Which[
    Abs[d . wHat - 1] < 0.0001,  IdentityMatrix[4],
    Abs[d . wHat + 1] < 0.0001, -IdentityMatrix[4],
    True,
      ang = ArcCos[Clip[d . wHat, {-1, 1}]];
      RotationMatrix[ang, {d, wHat}]
  ]
];

(* 隣接リスト（全セル，1回だけ計算） *)
Print["Building adjacency list ..."];
t0 = AbsoluteTime[];
adjFull = Table[
  Select[Range[nc], j |-> j != i &&
    Length[Intersection[cells[[i]], cells[[j]]]] >= 1],
  {i, nc}
];
Print["  Done: ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* RZ ペーリング（3D-face-centroid-up，単一ペア） *)
peel1[c1_, c2_, rotDir_] := Module[
  {localVers, R, cC, c1Vec, c2Vec,
   adj, order, last, newCell, cans, rCans,
   wvals, maxw, gvals, tied, geoScore},

  localVers = # - Mean[vers] & /@ vers;
  R = faceUpRot[rotDir];
  localVers = R . # & /@ localVers;
  cC = Mean[localVers[[ cellVerts[#] ]]] & /@ Range[nc];

  c1Vec = cC[[c1]];
  c2Vec = cC[[c2]];

  adj = Complement[#, {c1}] & /@ adjFull;
  order = {c1, c2};
  adj = Complement[#, {c2}] & /@ adj;
  last = c2;

  geoScore[k_] := If[Length[order] == 2,
    c2Vec[[1]] * cC[[k, 2]] - c2Vec[[2]] * cC[[k, 1]],
    Quiet[Det[{c1Vec, c2Vec, cC[[last]], cC[[k]]}]]
  ];

  While[Length[Flatten[adj]] != 0 && Length[adj[[last]]] != 0,
    cans = If[Length[order] == 2,
      adj[[last]],
      Select[adj[[last]],
        j |-> Quiet[Det[{c1Vec, c2Vec, cC[[last]], cC[[j]]}]] >= -eps]
    ];
    newCell = If[Length[cans] >= 1,
      wvals = cC[[cans, 4]]; maxw = Max[wvals];
      tied  = Select[cans, cC[[#, 4]] >= maxw - 10^-10 &];
      If[Length[tied] == 1, tied[[1]],
        gvals = geoScore[#] & /@ tied;
        tied[[ First[Ordering[-gvals, 1]] ]]],
      rCans = adj[[last]];
      wvals = cC[[rCans, 4]]; maxw = Min[wvals];
      tied  = Select[rCans, cC[[#, 4]] <= maxw + 10^-10 &];
      If[Length[tied] == 1, tied[[1]],
        gvals = Abs[geoScore[#]] & /@ tied;
        tied[[ First[Ordering[-gvals, 1]] ]]]
    ];
    last = newCell;
    AppendTo[order, newCell];
    adj = Complement[#, {newCell}] & /@ adj;
  ];

  Length[order]   (* 配置できたセル数 *)
];

(* ── まず C1=1 の4ペアでタイミング計測 ── *)
Print["\nTiming test: C1=1, all C2 candidates ..."];
c2s1 = adjFull[[1]];
Print["  C2 candidates for C1=1: ", c2s1];
t0 = AbsoluteTime[];
Do[
  c2 = c2s1[[k]];
  sharedFI = Intersection[cells[[1]], cells[[c2]]][[1]];
  sharedFC = Mean[vers[[ faces[[sharedFI]] ]]];
  nPlaced  = peel1[1, c2, sharedFC];
  Print["  C1=1 C2=", c2, ": ", nPlaced, "/", nc,
        " (", If[nPlaced==nc,"SUCCESS","FAIL"], ")"],
  {k, Length[c2s1]}
];
tPer = (AbsoluteTime[] - t0) / Length[c2s1];
Print["  Per-pair time: ", Round[tPer, 0.1], " s"];
Print["  Estimated total (2400 pairs): ", Round[2400 * tPer / 60, 1], " min"];

(* ── 全 2400 ペア ── *)
Print["\nRunning all ", nc, " * 4 = 2400 pairs (3D-face-centroid-up + RZ) ..."];
t0 = AbsoluteTime[];
nOK = 0; nFail = 0;
stepDist = <||>;   (* 終了ステップ → 件数 *)

Do[
  Do[
    c2 = adjFull[[c1, k]];
    sharedFI = Intersection[cells[[c1]], cells[[c2]]][[1]];
    sharedFC = Mean[vers[[ faces[[sharedFI]] ]]];
    nPlaced  = peel1[c1, c2, sharedFC];
    If[nPlaced == nc,
      nOK++,
      nFail++;
      stepDist[nPlaced] =
        Lookup[stepDist, nPlaced, 0] + 1
    ],
    {k, Length[adjFull[[c1]]]}
  ];
  If[Mod[c1, 60] == 0 || c1 == nc,
    Print["  C1=", c1, "/", nc,
          "  OK=", nOK, "  Fail=", nFail,
          "  (", Round[AbsoluteTime[] - t0, 1], " s)"]
  ],
  {c1, nc}
];

tTotal = Round[AbsoluteTime[] - t0, 1];
Print["\n=== Results ==="];
Print["Total pairs: ", nOK + nFail];
Print["SUCCESS: ", nOK, " (", NumberForm[N[nOK/(nOK+nFail)*100], {4,1}], "%)"];
Print["FAIL:    ", nFail];
Print["Time:    ", tTotal, " s"];

If[nFail > 0,
  Print["\nFail step distribution:"];
  sorted = SortBy[Normal[stepDist], First];
  Do[
    Print["  step ", s[[1]], ": ", s[[2]], " pairs (",
          NumberForm[N[s[[2]]/(nOK+nFail)*100], {4,1}], "%)"],
    {s, sorted}
  ]
];

Print["Done."];
