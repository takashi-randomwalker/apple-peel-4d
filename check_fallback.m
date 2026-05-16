(* check_fallback.m
   目的:
     1) 4D ペーリングで fallback が発火するか計測 (peeling4DfLog)
     2) fallback なし版との成功率比較    (peeling4DfNoFB)
     3) 3D 正十二面体 R3 の正確な成功数を取得
*)

SetDirectory["/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D"];

(* ============================================================
   Part 1: fallback カウント付き 4D 関数
   ============================================================ *)

peeling4DfLog[origVers_, faces_, cells_, top_] :=
  Module[
    {wHat, cc1, localVers, cCenters,
     neighbors, neighborsInit,
     order, cansInit, chosenFrom1,
     lastChosen, newChoice, cans,
     angle4D, rCans, cellVerts, fbCount},

    wHat = {0, 0, 0, 1};
    cellVerts[i_] := Union @@ faces[[ cells[[i]] ]];
    localVers = (# - Mean[origVers]) & /@ origVers;
    cc1 = Mean[ localVers[[ cellVerts[top] ]] ];
    localVers = Which[
      Abs[(cc1/Norm[cc1]) . wHat - 1] < 0.0001, localVers,
      Abs[(cc1/Norm[cc1]) . wHat + 1] < 0.0001, -localVers,
      True,
        angle4D = ArcCos[Clip[cc1 . wHat / Norm[cc1], {-1,1}]];
        RotationMatrix[angle4D, {cc1/Norm[cc1], wHat}] . # & /@ localVers
    ];
    cCenters = Mean[ localVers[[ cellVerts[#] ]] ] & /@ Range[Length[cells]];
    neighbors = Table[
      Select[Range[Length[cells]],
        Function[j, j != i && Length[Intersection[cells[[i]], cells[[j]]]] >= 1]],
      {i, Length[cells]}
    ];
    neighbors = neighborsInit = Complement[#, {top}] & /@ neighbors;
    cansInit = neighbors[[top]];

    Table[
      order       = {top};
      neighbors   = neighborsInit;
      chosenFrom1 = cansInit[[i]];
      AppendTo[order, chosenFrom1];
      neighbors   = Complement[#, {chosenFrom1}] & /@ neighbors;
      lastChosen  = chosenFrom1;
      newChoice   = lastChosen;
      fbCount     = 0;

      While[
        Length[Flatten[neighbors]] != 0 &&
        Length[neighbors[[newChoice]]] != 0,

        cans = Select[
          neighbors[[lastChosen]],
          Function[j,
            With[{c = cCenters[[lastChosen]]},
              c[[2]] * cCenters[[j]][[1]] -
              c[[1]] * cCenters[[j]][[2]] >= 0]]
        ];

        newChoice = Which[
          Length[cans] > 1,
            cans[[ First@First@Position[cCenters[[cans,4]], Max[cCenters[[cans,4]]]] ]],
          Length[cans] == 1,
            cans[[1]],
          True,
            fbCount++;
            rCans = neighbors[[lastChosen]];
            rCans[[ First@First@Position[cCenters[[rCans,4]], Min[cCenters[[rCans,4]]]] ]]
        ];

        lastChosen = newChoice;
        AppendTo[order, newChoice];
        neighbors = Complement[#, {newChoice}] & /@ neighbors;
      ];

      {Length[order] == Length[cells], fbCount}
      , {i, Length[cansInit]}
    ]
  ];

(* ============================================================
   Part 1b: fallback なし版 (左候補なし → 即 Break)
   ============================================================ *)

peeling4DfNoFB[origVers_, faces_, cells_, top_] :=
  Module[
    {wHat, cc1, localVers, cCenters,
     neighbors, neighborsInit,
     order, cansInit, chosenFrom1,
     lastChosen, newChoice, cans,
     angle4D, cellVerts},

    wHat = {0, 0, 0, 1};
    cellVerts[i_] := Union @@ faces[[ cells[[i]] ]];
    localVers = (# - Mean[origVers]) & /@ origVers;
    cc1 = Mean[ localVers[[ cellVerts[top] ]] ];
    localVers = Which[
      Abs[(cc1/Norm[cc1]) . wHat - 1] < 0.0001, localVers,
      Abs[(cc1/Norm[cc1]) . wHat + 1] < 0.0001, -localVers,
      True,
        angle4D = ArcCos[Clip[cc1 . wHat / Norm[cc1], {-1,1}]];
        RotationMatrix[angle4D, {cc1/Norm[cc1], wHat}] . # & /@ localVers
    ];
    cCenters = Mean[ localVers[[ cellVerts[#] ]] ] & /@ Range[Length[cells]];
    neighbors = Table[
      Select[Range[Length[cells]],
        Function[j, j != i && Length[Intersection[cells[[i]], cells[[j]]]] >= 1]],
      {i, Length[cells]}
    ];
    neighbors = neighborsInit = Complement[#, {top}] & /@ neighbors;
    cansInit = neighbors[[top]];

    Table[
      order       = {top};
      neighbors   = neighborsInit;
      chosenFrom1 = cansInit[[i]];
      AppendTo[order, chosenFrom1];
      neighbors   = Complement[#, {chosenFrom1}] & /@ neighbors;
      lastChosen  = chosenFrom1;
      newChoice   = lastChosen;

      While[
        Length[Flatten[neighbors]] != 0 &&
        Length[neighbors[[newChoice]]] != 0,

        cans = Select[
          neighbors[[lastChosen]],
          Function[j,
            With[{c = cCenters[[lastChosen]]},
              c[[2]] * cCenters[[j]][[1]] -
              c[[1]] * cCenters[[j]][[2]] >= 0]]
        ];

        If[Length[cans] == 0, Break[]];   (* fallback なし: 即終了 *)

        newChoice = If[Length[cans] > 1,
          cans[[ First@First@Position[cCenters[[cans,4]], Max[cCenters[[cans,4]]]] ]],
          cans[[1]]
        ];

        lastChosen = newChoice;
        AppendTo[order, newChoice];
        neighbors = Complement[#, {newChoice}] & /@ neighbors;
      ];

      Length[order] == Length[cells]
      , {i, Length[cansInit]}
    ]
  ];

(* ============================================================
   Part 2: 4D 計算 (5/8/16/24-cell)
   ============================================================ *)

loadPoly[n_] := Module[{raw},
  raw = Get["_4DData/f" <> ToString[n] <> ".m"];
  <|"vers" -> N[raw[[1]]], "faces" -> raw[[4]], "cells" -> raw[[6]]|>
];

Print[""];
Print["==== 4D fallback check & no-fallback comparison ===="];
Print[StringRepeat["-", 76]];
Print[PaddedForm["Polytope",   12],
      PaddedForm["Pairs",       7],
      PaddedForm["WithFB",      9],
      PaddedForm["NoFB",        7],
      PaddedForm["Diff",        7],
      PaddedForm["FB-pairs",   10],
      PaddedForm["FB-total",   10]];
Print[StringRepeat["-", 76]];

Do[
  poly = loadPoly[n];
  nCells = Length[poly["cells"]];

  resLog = Flatten[
    Table[peeling4DfLog[poly["vers"], poly["faces"], poly["cells"], top],
          {top, nCells}], 1];
  resNoFB = Flatten[
    Table[peeling4DfNoFB[poly["vers"], poly["faces"], poly["cells"], top],
          {top, nCells}], 1];

  sucWithFB = Count[resLog,  {True, _}];
  sucNoFB   = Count[resNoFB, True];
  fbPairs   = Count[resLog,  {_, fb_} /; fb > 0];
  fbTotal   = Total[resLog[[All, 2]]];

  Print[PaddedForm[ToString[n]<>"-cell", 12],
        PaddedForm[Length[resLog],  7],
        PaddedForm[sucWithFB,       9],
        PaddedForm[sucNoFB,         7],
        PaddedForm[sucWithFB - sucNoFB, 7],
        PaddedForm[fbPairs,        10],
        PaddedForm[fbTotal,        10]];
  , {n, {5, 8, 16, 24}}
];
Print[StringRepeat["-", 76]];
Print["WithFB  : fallback あり版の成功数 (= 現行 peeling4Df.m)"];
Print["NoFB    : fallback なし版の成功数 (左候補なし → 即終了)"];
Print["Diff    : WithFB - NoFB  (fallback による上乗せ分)"];
Print["FB-pairs: fallback が1回以上発火した (C1,C2) ペア数"];
Print["FB-total: 全ペアの fallback 発火回数合計"];

(* ============================================================
   Part 3: 3D 正十二面体 R3
   ============================================================ *)

Print[""];
Print["==== 3D Dodecahedron R3 (maxz) ===="];
Get["peeling3DLoxo.m"];
dod = p3lFromBuiltin["Dodecahedron", "maxz"];
Print["complete (total pairs that succeeded) : ", dod["complete"]];
Print["unique (distinct orderings)           : ", dod["unique"]];
Print["total (C1,C2) pairs evaluated         : ", dod["total"]];
sucPct = N[100 * dod["complete"] / dod["total"], 4];
Print["success rate                          : ", sucPct, " %"];
