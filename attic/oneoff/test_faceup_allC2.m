(* ================================================================
   test_faceup_allC2.m
   C1=3 の全 C2 候補（10通り）に対して
     A: cell-centroid-up  (現行)
     B: 3D-face-centroid-up  (共有面の重心を +w に揃える)
   の2方式で RZ ペーリング → SAT 重なり判定し結果を比較する．
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];
Get[scriptDir <> "unfold3DExport.m"];

dataDir = scriptDir <> "_4DData/";
raw   = Get[dataDir <> "f120.m"];
vers  = N[raw[[1]]]; faces = raw[[4]]; cells = raw[[6]];
nc    = Length[cells];

fixedC1 = 3;
(* C1=3 の成功 C2 一覧（analyze_120cell_c2types.m で確認済み） *)
successC2list = {1, 2, 4, 10, 15, 17, 18, 97, 100, 115};

cellVerts[ci_] := Union @@ faces[[ cells[[ci]] ]];
wHat = {0., 0., 0., 1.};

faceUpRot[dir_] := Module[{d = Normalize[dir - Mean[vers]], ang},
  Which[
    Abs[d . wHat - 1] < 0.0001,  IdentityMatrix[4],
    Abs[d . wHat + 1] < 0.0001, -IdentityMatrix[4],
    True,
      ang = ArcCos[Clip[d . wHat, {-1, 1}]];
      RotationMatrix[ang, {d, wHat}]
  ]
];

(* 単一ペアのペーリング（RZ）．rotDir = face-up の方向ベクトル（元座標） *)
peel1Pair[c1_, c2_, rotDir_] := Module[
  {eps, localVers, R, cC, c1Vec, c2Vec,
   adjList, order, last, newCell, cans, rCans,
   wvals, maxw, gvals, tied, geoScore},

  eps = 10^-10;
  localVers = # - Mean[vers] & /@ vers;
  R = faceUpRot[rotDir];
  localVers = R . # & /@ localVers;
  cC = Mean[localVers[[ cellVerts[#] ]]] & /@ Range[nc];

  c1Vec = cC[[c1]];
  c2Vec = cC[[c2]];

  adjList = Table[
    Select[Range[nc], j |-> j != i &&
      Length[Intersection[cells[[i]], cells[[j]]]] >= 1],
    {i, nc}
  ];
  adjList = Complement[#, {c1}] & /@ adjList;

  order = {c1, c2};
  adjList = Complement[#, {c2}] & /@ adjList;
  last = c2;

  geoScore[k_] := If[Length[order] == 2,
    c2Vec[[1]] * cC[[k, 2]] - c2Vec[[2]] * cC[[k, 1]],
    Quiet[Det[{c1Vec, c2Vec, cC[[last]], cC[[k]]}]]
  ];

  While[Length[Flatten[adjList]] != 0 && Length[adjList[[last]]] != 0,
    cans = If[Length[order] == 2,
      adjList[[last]],
      Select[adjList[[last]],
        j |-> Quiet[Det[{c1Vec, c2Vec, cC[[last]], cC[[j]]}]] >= -eps]
    ];
    newCell = If[Length[cans] >= 1,
      wvals = cC[[cans, 4]]; maxw = Max[wvals];
      tied  = Select[cans, cC[[#, 4]] >= maxw - 10^-10 &];
      If[Length[tied] == 1, tied[[1]],
        gvals = geoScore[#] & /@ tied;
        tied[[ First[Ordering[-gvals, 1]] ]]],
      rCans = adjList[[last]];
      wvals = cC[[rCans, 4]]; maxw = Min[wvals];
      tied  = Select[rCans, cC[[#, 4]] <= maxw + 10^-10 &];
      If[Length[tied] == 1, tied[[1]],
        gvals = Abs[geoScore[#]] & /@ tied;
        tied[[ First[Ordering[-gvals, 1]] ]]]
    ];
    last = newCell;
    AppendTo[order, newCell];
    adjList = Complement[#, {newCell}] & /@ adjList;
  ];

  {order, Length[order] == nc}
];

(* cellLocal3Ds を事前計算 *)
Print["Pre-computing cellLocal3Ds ..."];
t0 = AbsoluteTime[];
cellLocal3Ds = Table[unfoldCellLocal3D[ci, vers, cells, faces], {ci, nc}];
Print["  Done: ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* C1 の cell-centroid-up 方向（現行，C2 によらず共通） *)
c1CentroidDir = Mean[vers[[ cellVerts[fixedC1] ]]];

(* ── 全 C2 ループ ─────────────────────────────────────── *)
Print["\n", StringRepeat["-", 70]];
Print[StringPadLeft["C2",4],"  ",
      StringPadLeft["A-ok",5],"  ",StringPadLeft["A-ov",5],"  ",
      StringPadLeft["B-ok",5],"  ",StringPadLeft["B-ov",5],"  ",
      "  A-order=B?  A-order-len  B-order-len"];
Print[StringRepeat["-", 70]];

totA = 0; totB = 0;
Do[
  c2 = successC2list[[idx]];
  t1 = AbsoluteTime[];

  (* A: cell-centroid-up *)
  {ordA, okA} = peel1Pair[fixedC1, c2, c1CentroidDir];
  If[okA,
    embsA = unfoldTo3DFast[cells, faces, ordA, cellLocal3Ds];
    ovsA  = Length[checkOverlaps[ordA, cells, faces, embsA]];
    , ovsA = -1];

  (* B: 3D-face-centroid-up（共有面の重心） *)
  sharedFI = Intersection[cells[[fixedC1]], cells[[c2]]][[1]];
  sharedFC = Mean[vers[[ faces[[sharedFI]] ]]];
  {ordB, okB} = peel1Pair[fixedC1, c2, sharedFC];
  If[okB,
    embsB = unfoldTo3DFast[cells, faces, ordB, cellLocal3Ds];
    ovsB  = Length[checkOverlaps[ordB, cells, faces, embsB]];
    , ovsB = -1];

  totA += If[ovsA >= 0, ovsA, 0];
  totB += If[ovsB >= 0, ovsB, 0];

  sameQ = If[okA && okB, ordA === ordB, "-"];
  Print[StringPadLeft[ToString[c2], 4], "  ",
        StringPadLeft[If[okA,"OK","FAIL"], 5], "  ",
        StringPadLeft[ToString[ovsA],     5], "  ",
        StringPadLeft[If[okB,"OK","FAIL"], 5], "  ",
        StringPadLeft[ToString[ovsB],     5], "  ",
        "  ", sameQ, "  ",
        Length[ordA], "  ", Length[ordB],
        "  (", Round[AbsoluteTime[] - t1, 1], " s)"],
  {idx, Length[successC2list]}
];

Print[StringRepeat["-", 70]];
Print["Total overlaps:  A = ", totA, "  B = ", totB,
      "  diff = ", totA - totB,
      If[totB < totA, "  (B wins by " <> ToString[totA-totB] <> ")",
         If[totA < totB, "  (A wins)", "  (tie)"]]];
Print["Done."];
