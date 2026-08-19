(* ================================================================
   analyze_600cell_stuck_example.m
   ================================================================
   600 胞体で C1=1 の全 4 つの C2 候補について peeling4Df4 (RZ) を
   実行し，各実行の終了時点における詰まりの構造を抽出する．

   各 (C1, C2) で：
     - 終了ステップ k (= 配置済みセル数)
     - 終了時の lastChosen セル id
     - lastChosen の 4 隣接セル（全グラフ上）
     - そのうち visited / unvisited の内訳
     - lastChosen の重心の w 座標
     - 直前の数ステップで visited を増やしていった経過

   論文の「Worked example: getting stuck」用の事実関係を抽出する．
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];

dataDir = scriptDir <> "_4DData/";
raw   = Get[dataDir <> "f600.m"];
vers  = raw[[1]];
faces = raw[[4]];
cells = raw[[6]];
nc    = Length[cells];
Print["600-cell loaded: ", nc, " cells."];

eps  = 10^-10;
wHat = {0, 0, 0, 1};

cellVerts[i_] := Union @@ faces[[ cells[[i]] ]];

(* Face-up alignment for C1=1 *)
top = 1;
localVers0 = (# - Mean[N[vers]]) & /@ N[vers];
cc1 = Mean[ localVers0[[ cellVerts[top] ]] ];
localVers = Which[
  Abs[(cc1 / Norm[cc1]) . wHat - 1] < 0.0001, localVers0,
  Abs[(cc1 / Norm[cc1]) . wHat + 1] < 0.0001, -localVers0,
  True,
    angle4D = ArcCos[Clip[cc1 . wHat / Norm[cc1], {-1, 1}]];
    RotationMatrix[angle4D, {cc1 / Norm[cc1], wHat}] . # & /@ localVers0
];
cCenters = Mean[ localVers[[ cellVerts[#] ]] ] & /@ Range[nc];
c1Vec    = cCenters[[top]];

Print["c1 (after Face-up) = ", N[c1Vec, 4]];
Print["c1 . wHat / |c1|  = ", N[(c1Vec.wHat) / Norm[c1Vec], 4]];

(* Adjacency list (all neighbors, ignoring visit state) *)
adjacencyAll = Table[
  Select[Range[nc],
    Function[j, j != i && Length[Intersection[cells[[i]], cells[[j]]]] >= 1]],
  {i, nc}
];

(* C2 candidates *)
cansInit = Complement[adjacencyAll[[top]], {top}];
Print["C2 candidates for C1=1: ", cansInit];

(* Instrumented peeling: runs RZ as in peeling4Df4.m, but also logs the
   final state at termination (no full trace to keep output small). *)
runOnePair[chosenFrom1_] := Module[
  {neighbors, order, lastChosen, newChoice, c2Vec, stuck,
   pending, cans, geoScore, gvals, wvals, maxw, tied, rCans,
   minw, terminalState},

  neighbors  = Complement[#, {top, chosenFrom1}] & /@ adjacencyAll;
  order      = {top, chosenFrom1};
  lastChosen = chosenFrom1;
  newChoice  = lastChosen;
  c2Vec      = cCenters[[chosenFrom1]];
  stuck      = False;

  While[
    !stuck &&
    Length[Flatten[neighbors]] != 0 &&
    Length[neighbors[[newChoice]]] != 0,

    pending = neighbors[[lastChosen]];
    cans = If[
      Length[order] == 2,
      pending,
      Select[pending,
        Function[j,
          Quiet[Det[{c1Vec, c2Vec, cCenters[[lastChosen]], cCenters[[j]]}]] >= -eps
        ]
      ]
    ];

    geoScore = If[Length[order] == 2,
      Function[j, c2Vec[[1]]*cCenters[[j,2]] - c2Vec[[2]]*cCenters[[j,1]]],
      Function[j, Quiet[Det[{c1Vec, c2Vec, cCenters[[lastChosen]], cCenters[[j]]}]]]
    ];

    newChoice = Which[
      Length[cans] >= 1,
        wvals = cCenters[[cans, 4]];
        maxw  = Max[wvals];
        tied  = Select[cans, cCenters[[#,4]] >= maxw - 10^-10 &];
        If[Length[tied] == 1,
          tied[[1]],
          gvals = geoScore[#] & /@ tied;
          If[Max[Abs[gvals]] < 10^-10,
            gvals = (c2Vec[[1]]*cCenters[[#,2]] - c2Vec[[2]]*cCenters[[#,1]]) & /@ tied
          ];
          tied[[ First[Ordering[-gvals, 1]] ]]
        ],
      True,
        rCans = pending;
        wvals = cCenters[[rCans, 4]];
        minw  = Min[wvals];
        tied  = Select[rCans, cCenters[[#,4]] <= minw + 10^-10 &];
        If[Length[tied] == 1,
          tied[[1]],
          gvals = Abs[geoScore[#]] & /@ tied;
          If[Max[gvals] < 10^-10,
            gvals = Abs[c2Vec[[1]]*cCenters[[#,2]] - c2Vec[[2]]*cCenters[[#,1]]] & /@ tied
          ];
          tied[[ First[Ordering[-gvals, 1]] ]]
        ]
    ];

    lastChosen = newChoice;
    AppendTo[order, newChoice];
    neighbors = Complement[#, {newChoice}] & /@ neighbors;
  ];

  (* Capture terminal state for diagnosis *)
  terminalState = <|
    "c2"             -> chosenFrom1,
    "orderLength"    -> Length[order],
    "success"        -> (Length[order] == nc),
    "lastCell"       -> lastChosen,
    "lastCell_w"     -> cCenters[[lastChosen, 4]],
    "lastCellAdjAll" -> adjacencyAll[[lastChosen]],
    "unvisited"      -> neighbors[[lastChosen]],
    "visitedFromAdj" -> Complement[adjacencyAll[[lastChosen]], neighbors[[lastChosen]]],
    "stepInOrder"    -> Position[order, lastChosen][[1, 1]]
  |>;

  terminalState
];

Print[""];
Print[StringRepeat["=", 70]];
Print["  Terminal state for each C2 (C1=1, rule=RZ, with fallback)"];
Print[StringRepeat["=", 70]];

results = runOnePair /@ cansInit;

Do[
  r = results[[i]];
  Print[""];
  Print["----- C2 = ", r["c2"], " -----"];
  Print["  order length      = ", r["orderLength"], "  (success: ", r["success"], ")"];
  Print["  last cell         = ", r["lastCell"], " (position ", r["stepInOrder"], " in order)"];
  Print["  w(last cell)      = ", NumberForm[r["lastCell_w"], {6, 4}]];
  Print["  all neighbors     = ", r["lastCellAdjAll"], "  (degree = ", Length[r["lastCellAdjAll"]], ")"];
  Print["  visited (in order)= ", r["visitedFromAdj"]];
  Print["  unvisited         = ", r["unvisited"], "  (size = ", Length[r["unvisited"]], ")"];
  , {i, Length[results]}
];

(* Pick one stuck example and dump full per-neighbor details:
   For each of the 4 neighbors of the last cell, show:
   - cell id
   - centroid w
   - whether visited (and its order index)
   - the Det value Det[{c1, c2, c_last, c_j}]
*)

Print[""];
Print[StringRepeat["=", 70]];
Print["  Detailed dump at the stuck step (first failing C2)"];
Print[StringRepeat["=", 70]];

failed = Select[results, !#["success"] &];
If[Length[failed] == 0,
  Print["  (no failure for C1=1 - unexpected)"];,
  example = failed[[1]];
  Print["  C2 = ", example["c2"], ", stuck at step k = ", example["orderLength"]];
  Print["  last cell        = ", example["lastCell"]];
  Print["  w(last)          = ", NumberForm[example["lastCell_w"], {6, 4}]];
  Print[""];
  (* re-run to get c2Vec for Det values *)
  c2v = cCenters[[example["c2"]]];
  cl  = example["lastCell"];
  Print["  Per-neighbor details (all 4 neighbors of last cell):"];
  Print["  ", StringPadRight["nbr id", 8],
        StringPadRight["w", 12],
        StringPadRight["Det[c1,c2,cL,cj]", 22],
        StringPadRight["visited?", 11],
        "order pos"];
  Print["  ", StringRepeat["-", 60]];
  Do[
    cw = cCenters[[nb, 4]];
    detVal = Quiet[Det[{c1Vec, c2v, cCenters[[cl]], cCenters[[nb]]}]];
    isVisited = MemberQ[example["visitedFromAdj"], nb];
    Print["  ",
      StringPadRight[ToString[nb], 8],
      StringPadRight[ToString[NumberForm[cw, {6, 4}]], 12],
      StringPadRight[ToString[NumberForm[detVal, {6, 4}]], 22],
      StringPadRight[If[isVisited, "yes", "no"], 11],
      If[isVisited, "(in order)", "-"]
    ],
    {nb, example["lastCellAdjAll"]}
  ];
];

Print["\nDone."];
