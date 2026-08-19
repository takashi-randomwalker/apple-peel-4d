(* ================================================================
   test_600cell_vertexup.m
   ================================================================
   600胞体で vertex-up（C1 の各頂点を +w に揃える）を試し，
   cell-centroid-up との停止ステップを比較する．
   C1=1, 全4頂点 × 全4 C2候補 = 16 ケースを計算．
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];
dataDir = scriptDir <> "_4DData/";
Get[scriptDir <> "peeling4Df4.m"];

raw   = Get[dataDir <> "f600.m"];
vers  = N[raw[[1]]];
faces = raw[[4]];
cells = raw[[6]];
nc    = Length[cells];

cellVerts[i_] := Union @@ faces[[ cells[[i]] ]];
wHat = {0, 0, 0, 1};
eps  = 10^-10;

(* ---- peelFromUp: up ベクトルを指定してペーリング（C2 固定） ---- *)
peelFromUp[upVec_, top_, c2Idx_] := Module[
  {lv, cc1, angle4D, cCenters, c1Vec, c2Vec,
   neighbors, order, lastChosen, newChoice, cans, stuck,
   geoScore, wvals, gvals, maxw, tied, rCans, minw},

  lv = (# - Mean[vers]) & /@ vers;

  (* up ベクトルを +w に整列 *)
  cc1 = upVec - Mean[vers];   (* 重心移動後の up ベクトル *)
  lv = Which[
    Abs[(cc1/Norm[cc1]).wHat - 1] < 0.0001, lv,
    Abs[(cc1/Norm[cc1]).wHat + 1] < 0.0001, -lv,
    True,
      angle4D = ArcCos[Clip[cc1.wHat/Norm[cc1], {-1,1}]];
      RotationMatrix[angle4D, {cc1/Norm[cc1], wHat}].# & /@ lv
  ];

  cCenters = Mean[lv[[ cellVerts[#] ]]] & /@ Range[nc];
  c1Vec = cCenters[[top]];
  c2Vec = cCenters[[c2Idx]];

  neighbors = Table[
    Select[Range[nc], j |-> j != i &&
      Length[Intersection[cells[[i]], cells[[j]]]] >= 1],
    {i, nc}
  ];
  neighbors = Complement[#, {top}] & /@ neighbors;

  order      = {top, c2Idx};
  neighbors  = Complement[#, {c2Idx}] & /@ neighbors;
  lastChosen = c2Idx;
  newChoice  = c2Idx;
  stuck      = False;

  geoScore[kIdx_, j_] := If[Length[order] == 2,
    c2Vec[[1]]*cCenters[[j,2]] - c2Vec[[2]]*cCenters[[j,1]],
    Quiet[Det[{c1Vec, c2Vec, cCenters[[kIdx]], cCenters[[j]]}]]
  ];

  While[!stuck && Length[Flatten[neighbors]] != 0 &&
        Length[neighbors[[newChoice]]] != 0,

    cans = If[Length[order] == 2,
      neighbors[[lastChosen]],
      Select[neighbors[[lastChosen]],
        j |-> Quiet[Det[{c1Vec, c2Vec, cCenters[[lastChosen]], cCenters[[j]]}]] >= -eps]
    ];

    If[Length[cans] == 0,
      (* fallback: min-w *)
      rCans = neighbors[[lastChosen]];
      wvals = cCenters[[rCans, 4]];
      minw  = Min[wvals];
      tied  = Select[rCans, cCenters[[#,4]] <= minw + eps &];
      newChoice = If[Length[tied] == 1, tied[[1]],
        gvals = Abs[geoScore[lastChosen, #]] & /@ tied;
        If[Max[gvals] < eps,
          gvals = Abs[c2Vec[[1]]*cCenters[[#,2]] - c2Vec[[2]]*cCenters[[#,1]]] & /@ tied];
        tied[[ First[Ordering[-gvals,1]] ]]
      ],
      (* main: max-w *)
      wvals = cCenters[[cans, 4]];
      maxw  = Max[wvals];
      tied  = Select[cans, cCenters[[#,4]] >= maxw - eps &];
      newChoice = If[Length[tied] == 1, tied[[1]],
        gvals = geoScore[lastChosen, #] & /@ tied;
        If[Max[Abs[gvals]] < eps,
          gvals = c2Vec[[1]]*cCenters[[#,2]] - c2Vec[[2]]*cCenters[[#,1]] & /@ tied];
        tied[[ First[Ordering[-gvals,1]] ]]
      ]
    ];

    lastChosen = newChoice;
    AppendTo[order, newChoice];
    neighbors = Complement[#, {newChoice}] & /@ neighbors;
  ];

  {Length[order], Length[order] == nc}
];

(* ---- C1=1 の隣接セル ---- *)
top = 1;
lv0 = (# - Mean[vers]) & /@ vers;
neighborsC1 = Select[Range[nc],
  j |-> j != top && Length[Intersection[cells[[top]], cells[[j]]]] >= 1];
verts1 = cellVerts[top];   (* C1 の頂点インデックス *)

Print["C1=", top, "  vertices: ", verts1, "  C2 candidates: ", neighborsC1];

(* ---- cell-centroid-up（現行）の参照ケース ---- *)
Print["\n--- cell-centroid-up (reference) ---"];
Do[
  res = peeling4Df4[vers, faces, cells, top, "RZ", False];
  c2res = res[[Position[neighborsC1, c2][[1,1]]]];
  Print["  C2=", c2, "  steps=", Length[c2res[[2]]], "  success=", c2res[[3]]],
  {c2, neighborsC1}
];

(* ---- vertex-up: C1 の各頂点で試す ---- *)
Print["\n--- vertex-up: all 4 vertices of C1=", top, " ---"];
Print[StringPadLeft["VertexUp", 10], "  ",
      StringPadLeft["C2", 5], "  ",
      StringPadLeft["Steps", 6], "  ",
      StringPadLeft["Success", 8]];
Print[StringRepeat["-", 35]];

Do[
  vIdx = verts1[[vi]];
  upVec = vers[[vIdx]];   (* この頂点を +w に揃える *)
  Do[
    {nSteps, ok} = peelFromUp[upVec, top, c2];
    Print[StringPadLeft["v" <> ToString[vi] <> "(#" <> ToString[vIdx] <> ")", 10], "  ",
          StringPadLeft[ToString[c2], 5], "  ",
          StringPadLeft[ToString[nSteps], 6], "  ",
          StringPadLeft[If[ok, "OK", "STUCK"], 8]],
    {c2, neighborsC1}
  ],
  {vi, 4}
];

Print["\nDone."];
