(* ================================================================
   peeling4Df4.m
   ================================================================
   4次元正多胞体 Apple-Peel ペーリング -- 規則3 + c1-c2 平面グローバル参照

   左条件（4次元体積形式）:
     Det[{c_1, c_2, c_k, c_j}] >= -eps

   3D 版との対応:
     3D: Det[{c_1, c_k, c_j}]        >= -eps  (c_1 固定, 北極 = +z)
     4D: Det[{c_1, c_2, c_k, c_j}]   >= -eps  (c_1–c_2 平面固定, 北極 = +w)

   Face-up 後（c_1 = (0,0,0,w_1)）の代数展開:
     = -w_1 * Det_xyz[{c_2, c_k, c_j}]

   k=2 の特別扱い不要:
     c_k = c_2 のとき Det = 0 → 全候補通過

   等変性:
     A = FaceUp(C1') ∘ σ ∘ FaceUp(C1)^{-1} は SO(4) で det(A)=1
     → Det 条件は σ に関して等変（k=2 を除く）

   幾何スコア（geo-score）:
     k=2:  c2_x*cj_y - c2_y*cj_x  （xy 平面符号付き面積）
     k≥3:  Det[{c1, c2, c_k, c_j}] （4次元体積形式）

   選択規則（rule 引数で切り替え）:
     RZ（default）:
       主基準 = max w；二次 = max geo-score
       フォールバック = min w；二次 = max |geo-score|
     RS:
       主基準 = max geo-score；二次 = max w
       フォールバック = max |geo-score|；二次 = min w

   eps = 10^-10

   引数:
     peeling4Df4[origVers, faces, cells, top, rule:"RZ"]
     rule = "RZ" | "RS"
   戻り値:
     { {localVers, order, successBool}, ... }  C2 候補ごとに1要素
   ================================================================ *)

Clear[peeling4Df4];

peeling4Df4[origVers_, faces_, cells_, top_, rule_:"RZ", noFB_:False] :=
  Module[
    {eps, wHat, cc1, localVers, cCenters,
     neighbors, neighborsInit,
     order, cansInit, chosenFrom1,
     lastChosen, newChoice, cans, stuck,
     c1Vec, c2Vec, angle4D, cellVerts,
     geoScore, wvals, gvals, maxw, maxg, maxAbsG, minw, tied, rCans},

    eps  = 10^-10;
    wHat = {0, 0, 0, 1};

    cellVerts[i_] := Union @@ faces[[ cells[[i]] ]];

    (* Step 1: 重心を原点へ平行移動 *)
    localVers = (# - Mean[origVers]) & /@ origVers;

    (* Step 2: 出発セル C1 の重心 *)
    cc1 = Mean[ localVers[[ cellVerts[top] ]] ];

    (* Step 3: C1 を +w 方向に整列（Face-up） *)
    localVers = Which[
      Abs[(cc1 / Norm[cc1]) . wHat - 1] < 0.0001,
        localVers,
      Abs[(cc1 / Norm[cc1]) . wHat + 1] < 0.0001,
        -localVers,
      True,
        angle4D = ArcCos[Clip[cc1 . wHat / Norm[cc1], {-1, 1}]];
        RotationMatrix[angle4D, {cc1 / Norm[cc1], wHat}] . # & /@ localVers
    ];

    (* Step 4: 各セルの重心（Face-up 回転後） *)
    cCenters = Mean[ localVers[[ cellVerts[#] ]] ] & /@ Range[Length[cells]];

    (* C1 のグローバル参照ベクトル（全ペア共通） *)
    c1Vec = cCenters[[top]];

    (* Step 5: 隣接リスト *)
    neighbors = Table[
      Select[
        Range[Length[cells]],
        Function[j,
          j != i && Length[Intersection[cells[[i]], cells[[j]]]] >= 1]
      ],
      {i, Length[cells]}
    ];

    (* Step 6: C1 を隣接リストから除去 *)
    neighbors = neighborsInit = Complement[#, {top}] & /@ neighbors;

    (* Step 7: C2 候補 = C1 の隣接セル *)
    cansInit = neighbors[[top]];

    (* Step 8: C2 ごとのペーリングループ *)
    Table[
      order       = {top};
      neighbors   = neighborsInit;
      chosenFrom1 = cansInit[[i]];
      AppendTo[order, chosenFrom1];
      neighbors   = Complement[#, {chosenFrom1}] & /@ neighbors;
      lastChosen  = chosenFrom1;
      newChoice   = lastChosen;
      c2Vec       = cCenters[[chosenFrom1]];
      stuck       = False;

      While[
        !stuck &&
        Length[Flatten[neighbors]] != 0 &&
        Length[neighbors[[newChoice]]] != 0,

        (* ── 左半空間フィルタ ──────────────────────────────── *)
        cans = If[
          Length[order] == 2,
          neighbors[[lastChosen]],     (* k=2: Det=0 → 全候補通過 *)
          Select[
            neighbors[[lastChosen]],
            Function[j,
              Quiet[Det[{c1Vec, c2Vec,
                         cCenters[[lastChosen]], cCenters[[j]]}]] >= -eps
            ]
          ]
        ];

        (* ── 幾何スコア（geo-score）を計算する補助 ─────────── *)
        (* k=2: xy 平面符号付き面積；k≥3: Det *)
        geoScore = If[Length[order] == 2,
          Function[j,
            c2Vec[[1]] * cCenters[[j, 2]] - c2Vec[[2]] * cCenters[[j, 1]]
          ],
          Function[j,
            Quiet[Det[{c1Vec, c2Vec, cCenters[[lastChosen]], cCenters[[j]]}]]
          ]
        ];

        (* ── フォールバックなし時の早期終了 ────────────────── *)
        If[noFB && Length[cans] == 0, stuck = True; Break[]];

        (* ── 選択規則 ──────────────────────────────────────── *)
        newChoice = Which[

          (* 左候補あり *)
          Length[cans] >= 1,
          If[rule === "RZ",
            (* RZ: 主 = max w，二次 = max geo-score *)
            wvals = cCenters[[cans, 4]];
            maxw  = Max[wvals];
            tied  = Select[cans, cCenters[[#, 4]] >= maxw - 10^-10 &];
            If[Length[tied] == 1,
              tied[[1]],
              gvals = geoScore[#] & /@ tied;
              If[Max[Abs[gvals]] < 10^-10,
                gvals = (c2Vec[[1]]*cCenters[[#,2]] - c2Vec[[2]]*cCenters[[#,1]]) & /@ tied
              ];
              tied[[ First[Ordering[-gvals, 1]] ]]
            ],
            (* RS: 主 = max geo-score，二次 = max w *)
            gvals = geoScore[#] & /@ cans;
            maxg  = Max[gvals];
            tied  = Pick[cans, Map[# >= maxg - 10^-10 &, gvals]];
            If[Length[tied] == 1,
              tied[[1]],
              tied[[ First[Ordering[-cCenters[[tied, 4]], 1]] ]]
            ]
          ],

          (* フォールバック（左候補なし） *)
          True,
          rCans = neighbors[[lastChosen]];
          If[rule === "RZ",
            (* RZ fallback: 主 = min w，二次 = max |geo-score| *)
            wvals = cCenters[[rCans, 4]];
            minw  = Min[wvals];
            tied  = Select[rCans, cCenters[[#, 4]] <= minw + 10^-10 &];
            If[Length[tied] == 1,
              tied[[1]],
              gvals = Abs[geoScore[#]] & /@ tied;
              If[Max[gvals] < 10^-10,
                gvals = Abs[(c2Vec[[1]]*cCenters[[#,2]] - c2Vec[[2]]*cCenters[[#,1]])] & /@ tied
              ];
              tied[[ First[Ordering[-gvals, 1]] ]]
            ],
            (* RS fallback: 主 = max |geo-score|，二次 = min w *)
            gvals   = Abs[geoScore[#]] & /@ rCans;
            maxAbsG = Max[gvals];
            tied    = Pick[rCans, Map[# >= maxAbsG - 10^-10 &, gvals]];
            If[Length[tied] == 1,
              tied[[1]],
              tied[[ First[Ordering[cCenters[[tied, 4]], 1]] ]]
            ]
          ]
        ];

        lastChosen = newChoice;
        AppendTo[order, newChoice];
        neighbors = Complement[#, {newChoice}] & /@ neighbors;
      ];

      {localVers, order, !stuck && Length[order] == Length[cells]}
      , {i, Length[cansInit]}
    ]

  ];

Print["peeling4Df4.m loaded."];
Print["Left condition: Det[{c1,c2,c_k,c_j}] >= -eps  (4D volume form)"];
Print["  k=2: all pass (Det=0); geo-score = c2_x*cj_y - c2_y*cj_x"];
Print["  k>=3 Det=0 fallback: geo-score -> xy cross product (same as k=2)"];
Print["  RZ: max-w primary, max-geo secondary; fallback min-w, max|geo|"];
Print["  RS: max-geo primary, max-w secondary; fallback max|geo|, min-w"];
