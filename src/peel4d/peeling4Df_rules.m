(* ================================================================
   peeling4Df_rules.m
   ================================================================
   4次元正多胞体における Apple-Peel ペーリング（三規則対応版）

   規則:
     "maxphi" (R1) : 左候補から max φ (最も強い左旋回)
     "loxo"   (R2) : 左候補から min φ (loxodrome 近似)
     "maxw"   (R3) : 左候補から max w (剥き軸優先)  ← peeling4Df.m と同じ

   左半空間条件 (全規則共通):
     c_{k,y} * c'_x - c_{k,x} * c'_y >= 0   (xy 2成分外積)

   Darboux フレーム (4D 版, xy 射影):
     fHat : xy 平面への (c_k - c_{k-1}) の射影を正規化
     lHat : fHat を 90° 反時計回り = {-fHat[[2]], fHat[[1]]}
     φ    : c_j の xy 射影と fHat のなす角 (atan2 で符号付き)

   引数:
     origVers : 頂点座標 {{x,y,z,w}, ...}
     faces    : 面（頂点インデックスリスト）
     cells    : セル（面インデックスリスト）
     top      : 開始セル F1 のインデックス
     rule     : "maxphi" | "loxo" | "maxw"  (省略時 "maxw")

   戻り値: 各 F2 候補に対する {rotatedVers, order, successBool}
   ================================================================ *)


(* ----------------------------------------------------------------
   補助関数: xy 射影 Darboux フレーム
   ---------------------------------------------------------------- *)

(* xy 平面への前進・左ベクトルを返す (2次元) *)
p4rXYFrame[cPrev_, cCurr_] :=
  Module[{dxy, n},
    dxy = {cCurr[[1]] - cPrev[[1]], cCurr[[2]] - cPrev[[2]]};
    n   = Norm[dxy];
    If[n > 10^-10,
      {dxy / n, {-dxy[[2]], dxy[[1]]} / n},   (* fHat, lHat *)
      {{1., 0.}, {0., 1.}}                      (* フォールバック *)
    ]
  ];

(* xy 射影での角度 φ (符号付き, ラジアン)
   正 = 左, 負 = 右 *)
p4rAngleXY[fHat_, lHat_, cj_] :=
  Module[{proj, pn},
    proj = {cj[[1]], cj[[2]]};
    pn   = Norm[proj];
    If[pn < 10^-10, Return[0.]];
    proj = proj / pn;
    ArcTan[Clip[proj . fHat, {-1., 1.}],
           Clip[proj . lHat, {-1., 1.}]]
  ];


(* ----------------------------------------------------------------
   メイン関数
   ---------------------------------------------------------------- *)

peeling4DfRules[origVers_, faces_, cells_, top_, rule_:"maxw"] :=
  Module[
    {wHat, cc1, localVers, cCenters,
     neighbors, neighborsInit,
     order, cansInit, chosenFrom1,
     prevChosen, lastChosen, newChoice,
     cans, rCans, fHat, lHat, angles, absAngles,
     angle4D, cellVerts},

    wHat = {0, 0, 0, 1};

    (* cellVerts: セル i の頂点インデックス (面経由) *)
    cellVerts[i_] := Union @@ faces[[ cells[[i]] ]];

    (* Step 1: 重心を原点へ平行移動 *)
    localVers = (# - Mean[origVers]) & /@ origVers;

    (* Step 2: 開始セル C1 の重心 *)
    cc1 = Mean[ localVers[[ cellVerts[top] ]] ];

    (* Step 3: cc1 を +w 軸に整列 *)
    localVers = Which[
      Abs[(cc1 / Norm[cc1]) . wHat - 1] < 0.0001,
        localVers,
      Abs[(cc1 / Norm[cc1]) . wHat + 1] < 0.0001,
        -localVers,
      True,
        angle4D = ArcCos[Clip[cc1 . wHat / Norm[cc1], {-1, 1}]];
        RotationMatrix[angle4D, {cc1 / Norm[cc1], wHat}] . # & /@ localVers
    ];

    (* Step 4: 回転後のセル重心 *)
    cCenters = Mean[ localVers[[ cellVerts[#] ]] ] & /@ Range[Length[cells]];

    (* Step 5: 隣接リスト (面インデックスを 1 つ以上共有) *)
    neighbors = Table[
      Select[
        Range[Length[cells]],
        Function[j,
          j != i && Length[Intersection[cells[[i]], cells[[j]]]] >= 1]
      ],
      {i, Length[cells]}
    ];

    (* Step 6: F1 を除去 *)
    neighbors = neighborsInit = Complement[#, {top}] & /@ neighbors;

    (* Step 7: F2 候補 = F1 の隣接セル *)
    cansInit = neighbors[[top]];

    (* Step 8: 各 F2 についてペーリング *)
    Table[
      order       = {top};
      neighbors   = neighborsInit;
      chosenFrom1 = cansInit[[i]];
      AppendTo[order, chosenFrom1];
      neighbors   = Complement[#, {chosenFrom1}] & /@ neighbors;
      prevChosen  = top;
      lastChosen  = chosenFrom1;
      newChoice   = lastChosen;

      While[
        Length[Flatten[neighbors]] != 0 &&
        Length[neighbors[[newChoice]]] != 0,

        (* 左候補: xy 外積 >= 0 *)
        cans = Select[
          neighbors[[lastChosen]],
          Function[j,
            With[{c = cCenters[[lastChosen]]},
              c[[2]] * cCenters[[j]][[1]] -
              c[[1]] * cCenters[[j]][[2]] >= 0
            ]
          ]
        ];

        (* Darboux フレーム (xy 射影) *)
        {fHat, lHat} = p4rXYFrame[cCenters[[prevChosen]], cCenters[[lastChosen]]];

        newChoice = Which[

          (* R1: 左候補から max φ *)
          rule === "maxphi" && Length[cans] > 1,
            angles = p4rAngleXY[fHat, lHat, cCenters[[#]]] & /@ cans;
            cans[[ First[Ordering[-angles, 1]] ]],

          (* R2: 左候補から min φ *)
          rule === "loxo" && Length[cans] > 1,
            angles = p4rAngleXY[fHat, lHat, cCenters[[#]]] & /@ cans;
            cans[[ First[Ordering[angles, 1]] ]],

          (* R3: 左候補から max w *)
          rule === "maxw" && Length[cans] > 1,
            cans[[
              First[First[
                Position[cCenters[[cans, 4]], Max[cCenters[[cans, 4]]]]
              ]]
            ]],

          (* 左候補が 1 つ: 規則によらず選ぶ *)
          Length[cans] == 1,
            cans[[1]],

          (* フォールバック: 左候補なし *)
          True,
            rCans = neighbors[[lastChosen]];
            Which[
              rule === "maxphi" && Length[rCans] > 0,
                absAngles = Abs[p4rAngleXY[fHat, lHat, cCenters[[#]]]] & /@ rCans;
                rCans[[ First[Ordering[-absAngles, 1]] ]],

              rule === "loxo" && Length[rCans] > 0,
                absAngles = Abs[p4rAngleXY[fHat, lHat, cCenters[[#]]]] & /@ rCans;
                rCans[[ First[Ordering[absAngles, 1]] ]],

              rule === "maxw" && Length[rCans] > 0,
                rCans[[
                  First[First[
                    Position[cCenters[[rCans, 4]], Min[cCenters[[rCans, 4]]]]
                  ]]
                ]],

              True,
                Break[]; Null
            ]
        ];

        If[newChoice === Null, Break[]];
        prevChosen = lastChosen;
        lastChosen = newChoice;
        AppendTo[order, newChoice];
        neighbors = Complement[#, {newChoice}] & /@ neighbors;
      ];

      {localVers, order, Length[order] == Length[cells]}
      , {i, Length[cansInit]}
    ]

  ];


Print["peeling4Df_rules.m loaded."];
Print["Usage:"];
Print["  res = peeling4DfRules[vers, faces, cells, top, rule]"];
Print["  rule: \"maxphi\" (R1) | \"loxo\" (R2) | \"maxw\" (R3)"];
