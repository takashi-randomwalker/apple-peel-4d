(* ================================================================
   test_faceup_pair.m
   120胞体の単一ペア (C1, C2) に対して
     (A) cell-centroid-up  (現行)
     (B) 3D-face-centroid-up  (C1–C2 共有面の重心を +w に揃える)
   の2種類の face-up で RZ ペーリングを実行し，結果を比較する．
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
fixedC2 = 1;

cellVerts[ci_] := Union @@ faces[[ cells[[ci]] ]];
wHat = {0., 0., 0., 1.};

(* 共有面のインデックスと重心 *)
sharedFaceIdx = Intersection[cells[[fixedC1]], cells[[fixedC2]]];
Print["C1=", fixedC1, "  C2=", fixedC2,
      "  Shared face index: ", sharedFaceIdx];
sharedFC = Mean[vers[[ faces[[ sharedFaceIdx[[1]] ]] ]]];
Print["Shared face centroid: ", NumberForm[sharedFC, {6,4}]];

(* ── 共通サブルーチン ──────────────────────────────────── *)
faceUpRot[dir_] := Module[{d = Normalize[dir], ang},
  Which[
    Abs[d . wHat - 1] < 0.0001,  IdentityMatrix[4],
    Abs[d . wHat + 1] < 0.0001, -IdentityMatrix[4],
    True,
      ang = ArcCos[Clip[d . wHat, {-1, 1}]];
      RotationMatrix[ang, {d, wHat}]
  ]
];

peel1Pair[rotDir_, label_] := Module[
  {eps, localVers, R, cC, c1Vec, c2Vec,
   adjList, order, last, newCell, cans, rCans,
   wvals, maxw, gvals, maxg, tied, geoScore,
   steps},

  eps = 10^-10;

  (* 頂点を重心へ平行移動 → face-up 回転 *)
  localVers = # - Mean[vers] & /@ vers;
  R = faceUpRot[rotDir - Mean[vers]];
  localVers = R . # & /@ localVers;

  (* 全セル重心 *)
  cC = Mean[localVers[[ cellVerts[#] ]]] & /@ Range[nc];

  Print["\n=== ", label, " ==="];
  Print["  C1 centroid (face-up): ", NumberForm[cC[[fixedC1]], {6,4}]];
  Print["  C2 centroid (face-up): ", NumberForm[cC[[fixedC2]], {6,4}]];

  (* w 分布の確認：C1 の全隣接 C2 候補 *)
  c2cands = Select[Range[nc], # != fixedC1 &&
    Length[Intersection[cells[[fixedC1]], cells[[#]]]] >= 1 &];
  wC2 = Sort[cC[[c2cands, 4]]];
  Print["  C2 candidates w-range: [", Min[wC2]//N, ", ", Max[wC2]//N, "]  spread=",
        (Max[wC2]-Min[wC2])//N];
  Print["  Distinct w bands: ", Length[Union[Round[wC2, 0.001]]]];

  (* 参照ベクトル *)
  c1Vec = cC[[fixedC1]];
  c2Vec = cC[[fixedC2]];

  (* 隣接リスト（C1 を除外） *)
  adjList = Table[
    Select[Range[nc], j |-> j != i &&
      Length[Intersection[cells[[i]], cells[[j]]]] >= 1],
    {i, nc}
  ];
  adjList = Complement[#, {fixedC1}] & /@ adjList;

  (* ペーリング本体 *)
  order = {fixedC1, fixedC2};
  adjList = Complement[#, {fixedC2}] & /@ adjList;
  last = fixedC2;

  geoScore[k_] := If[Length[order] == 2,
    c2Vec[[1]] * cC[[k, 2]] - c2Vec[[2]] * cC[[k, 1]],
    Quiet[Det[{c1Vec, c2Vec, cC[[last]], cC[[k]]}]]
  ];

  While[
    Length[Flatten[adjList]] != 0 &&
    Length[adjList[[last]]] != 0,

    (* 左半空間フィルタ *)
    cans = If[Length[order] == 2,
      adjList[[last]],
      Select[adjList[[last]],
        j |-> Quiet[Det[{c1Vec, c2Vec, cC[[last]], cC[[j]]}]] >= -eps]
    ];

    (* RZ: max w → max geoScore *)
    newCell = If[Length[cans] >= 1,
      wvals = cC[[cans, 4]];
      maxw  = Max[wvals];
      tied  = Select[cans, cC[[#, 4]] >= maxw - 10^-10 &];
      If[Length[tied] == 1, tied[[1]],
        gvals = geoScore[#] & /@ tied;
        tied[[ First[Ordering[-gvals, 1]] ]]
      ],
      (* フォールバック: min w *)
      rCans = adjList[[last]];
      wvals = cC[[rCans, 4]];
      maxw  = Min[wvals];
      tied  = Select[rCans, cC[[#, 4]] <= maxw + 10^-10 &];
      If[Length[tied] == 1, tied[[1]],
        gvals = Abs[geoScore[#]] & /@ tied;
        tied[[ First[Ordering[-gvals, 1]] ]]
      ]
    ];

    last = newCell;
    AppendTo[order, newCell];
    adjList = Complement[#, {newCell}] & /@ adjList;
  ];

  success = Length[order] == nc;
  Print["  Result: ", If[success, "SUCCESS", "FAIL"],
        "  (", Length[order], "/", nc, " cells placed)"];
  {order, success}
];

(* ── (A) cell-centroid-up ────────────────────────────── *)
{ordA, okA} = peel1Pair[vers[[ cellVerts[fixedC1][[1]] ]] * 0 +
                         Mean[vers[[ cellVerts[fixedC1] ]]],
                         "A: cell-centroid-up"];

(* ── (B) 3D-face-centroid-up ─────────────────────────── *)
{ordB, okB} = peel1Pair[sharedFC,
                         "B: 3D-face-centroid-up (shared face)"];

Print["\n--- Comparison ---"];
Print["Same order? ", ordA === ordB];
If[okA && okB,
  Print["Both succeed. Overlap check (A):"];
  (* 重なり判定 *)
  cellLocal3Ds = Table[unfoldCellLocal3D[ci, vers, cells, faces], {ci, nc}];
  embsA = unfoldTo3DFast[cells, faces, ordA, cellLocal3Ds];
  ovsA  = checkOverlaps[ordA, cells, faces, embsA];
  Print["  Overlaps (A): ", Length[ovsA]];
  If[ordA =!= ordB,
    Print["Overlap check (B):"];
    embsB = unfoldTo3DFast[cells, faces, ordB, cellLocal3Ds];
    ovsB  = checkOverlaps[ordB, cells, faces, embsB];
    Print["  Overlaps (B): ", Length[ovsB]]
  ]
];

Print["Done."];
