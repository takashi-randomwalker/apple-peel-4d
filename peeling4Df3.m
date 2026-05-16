(* ================================================================
   peeling4Df3.m
   ================================================================
   4D Apple-Peel peeling -- Rule 3 (max-w) with proper left condition.

   差分 vs peeling4Df.m:
     左半空間条件を xy 2成分行列式から xyz 3成分行列式へ変更．
       旧: c_y * c'_x - c_x * c'_y >= 0
       新 (k>=2): Det[{c_{k-1}^xyz, c_k^xyz, c_j^xyz}] >= 0
       新 (k=1) : 旧条件を据え置き（直前セルが未定義のため）
     選択規則は変更なし: 左候補から max w，なければ min w（規則3）．

   引数・戻り値は peeling4Df と同一：
     peeling4Df3[origVers, faces, cells, top]
     -> { {localVers, order, successBool}, ... }  F2 候補ごとに1要素
   ================================================================ *)

peeling4Df3[origVers_, faces_, cells_, top_] :=
  Module[
    {wHat, cc1, localVers, cCenters,
     neighbors, neighborsInit,
     order, cansInit, chosenFrom1,
     lastChosen, prevChosen, newChoice, cans,
     angle4D, rCans, cellVerts},

    wHat = {0, 0, 0, 1};

    cellVerts[i_] := Union @@ faces[[ cells[[i]] ]];

    (* Step 1: 重心を原点に移動 *)
    localVers = (# - Mean[origVers]) & /@ origVers;

    (* Step 2: 出発セル C1 の重心 *)
    cc1 = Mean[ localVers[[ cellVerts[top] ]] ];

    (* Step 3: +w 方向へ回転 *)
    localVers = Which[
      Abs[(cc1 / Norm[cc1]) . wHat - 1] < 0.0001,
        localVers,
      Abs[(cc1 / Norm[cc1]) . wHat + 1] < 0.0001,
        -localVers,
      True,
        angle4D = ArcCos[Clip[cc1 . wHat / Norm[cc1], {-1, 1}]];
        RotationMatrix[angle4D, {cc1 / Norm[cc1], wHat}] . # & /@ localVers
    ];

    (* Step 4: 各セルの重心 (回転後) *)
    cCenters = Mean[ localVers[[ cellVerts[#] ]] ] & /@ Range[Length[cells]];

    (* Step 5: 隣接リスト (面インデックスを1つ以上共有) *)
    neighbors = Table[
      Select[
        Range[Length[cells]],
        Function[j,
          j != i && Length[Intersection[cells[[i]], cells[[j]]]] >= 1]
      ],
      {i, Length[cells]}
    ];

    (* Step 6: F1 を隣接リストから除去 *)
    neighbors = neighborsInit = Complement[#, {top}] & /@ neighbors;

    (* Step 7: F2 候補 = F1 の隣接セル *)
    cansInit = neighbors[[top]];

    (* Step 8: F2 ごとのペーリングループ *)
    Table[
      order       = {top};
      neighbors   = neighborsInit;
      chosenFrom1 = cansInit[[i]];
      AppendTo[order, chosenFrom1];
      neighbors   = Complement[#, {chosenFrom1}] & /@ neighbors;
      prevChosen  = top;           (* C1 が直前セル *)
      lastChosen  = chosenFrom1;   (* C2 が現在セル *)
      newChoice   = lastChosen;

      While[
        Length[Flatten[neighbors]] != 0 &&
        Length[neighbors[[newChoice]]] != 0,

        (* ---- 左半空間フィルタ ---- *)
        cans = Select[
          neighbors[[lastChosen]],
          Function[j,
            (* k=1 (最初のステップ) は xy 2成分条件で代替 *)
            If[Length[order] == 2,
              (* prevChosen = top, lastChosen = C2: xy条件 *)
              With[{c = cCenters[[lastChosen]]},
                c[[2]] * cCenters[[j]][[1]] -
                c[[1]] * cCenters[[j]][[2]] >= 0
              ],
              (* k>=2: xyz 3成分行列式条件 (定義 4.1) *)
              Det[{
                cCenters[[prevChosen]][[{1,2,3}]],
                cCenters[[lastChosen]][[{1,2,3}]],
                cCenters[[j]][[{1,2,3}]]
              }] >= 0
            ]
          ]
        ];

        (* ---- 規則3: 左候補から max w, なければ min w ---- *)
        newChoice = Which[
          Length[cans] > 1,
            cans[[ First[First[
              Position[cCenters[[cans, 4]], Max[cCenters[[cans, 4]]]]
            ]] ]],
          Length[cans] == 1,
            cans[[1]],
          True,
            rCans = neighbors[[lastChosen]];
            rCans[[ First[First[
              Position[cCenters[[rCans, 4]], Min[cCenters[[rCans, 4]]]]
            ]] ]]
        ];

        prevChosen = lastChosen;
        lastChosen = newChoice;
        AppendTo[order, newChoice];
        neighbors = Complement[#, {newChoice}] & /@ neighbors;
      ];

      {localVers, order, Length[order] == Length[cells]}
      , {i, Length[cansInit]}
    ]

  ];
