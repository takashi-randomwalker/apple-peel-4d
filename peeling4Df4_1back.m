(* ================================================================
   peeling4Df4_1back.m
   ================================================================
   1段バックトラック版 Apple-Peel ペーリング（4次元）

   左条件（peeling4Df4.m と同じ）:
     Det[{c_1, c_2, c_k, c_j}] >= -eps  (4次元体積形式，グローバル参照)
     k=2 は Det=0 → 全候補通過
     k=3 以降は c1–c2 平面に対する体積符号

   選択規則（peeling4Df4.m と同じ）:
     RZ (default): 主 = max w，二次 = max geo-score
                   fallback: min w，二次 = max |geo-score|
     RS:           主 = max geo-score，二次 = max w
                   fallback: max |geo-score|，二次 = min w

   バックトラック規則:
     詰まったとき（lastChosen の未試行隣接が空）:
       - justBT=True または Length[order]<=2 なら終了
       - そうでなければ: 直前1胞を取消し，tried に追加，状態復元
       - justBT = True にセット
     前進できたら justBT = False にリセット
     → 連続バックトラック不可・非連続バックトラック可

   インターフェース（peeling4Df4.m と同一）:
     peeling4Df4_1back[origVers, faces, cells, top, rule:"RZ"]
     戻り値: {{localVers, order, successBool}, ...}  C2 候補ごと

   eps = 10^-10
   ================================================================ *)

Clear[peeling4Df41back, p4f41backPair];

(* ── 単一 (top, f2) ペアの 1段バックトラック ─────────────────── *)
p4f41backPair[c1Vec_, c2Vec_, cCenters_, nc_,
              neighborsInit_, top_, f2_, rule_:"RZ"] :=
  Module[
    {eps, order, neighbors, lastChosen,
     tried, justBT, hasSaved, savedNbrs, savedLast,
     pendingRaw, triedHere, pending, cans, rCans,
     geoScore, wvals, gvals, maxw, maxg, maxAbsG, minw, tied, new, bad, newPos},

    eps = 10^-10;

    order      = {top, f2};
    neighbors  = Complement[#, {f2}] & /@ neighborsInit;
    lastChosen = f2;
    tried      = <||>;      (* pos -> {tried cells at that pos} *)
    justBT     = False;
    hasSaved   = False;
    savedNbrs  = {};
    savedLast  = 0;

    While[Length[order] < nc,

      (* 現在の lastChosen の未試行隣接 *)
      pendingRaw = neighbors[[lastChosen]];
      triedHere  = Lookup[tried, Length[order], {}];
      pending    = Complement[pendingRaw, triedHere];

      If[Length[pending] == 0,
        (* ─── 詰まり ─── *)
        If[justBT || Length[order] <= 2, Break[]];

        bad    = Last[order];
        order  = Most[order];
        If[hasSaved,
          neighbors  = savedNbrs;
          lastChosen = savedLast;
          hasSaved   = False
        ];
        newPos        = Length[order];
        tried[newPos] = Append[Lookup[tried, newPos, {}], bad];
        (* より後の位置のエントリを削除（古いパスのゴミを除去） *)
        tried  = Association @ KeySelect[tried, # <= newPos &];
        justBT = True,

        (* ─── 候補選択 ─── *)

        (* 左候補フィルタ（k=2 は Det=0 なので全通過） *)
        cans = If[
          Length[order] == 2,
          pending,
          Select[pending, Function[j,
            Quiet[Det[{c1Vec, c2Vec,
                       cCenters[[lastChosen]], cCenters[[j]]}]] >= -eps
          ]]
        ];

        (* フォールバック: 左候補なければ pending 全体 *)
        rCans = If[Length[cans] >= 1, cans, pending];

        (* geo-score 計算関数 *)
        geoScore = If[
          Length[order] == 2,
          Function[j,
            c2Vec[[1]] cCenters[[j, 2]] - c2Vec[[2]] cCenters[[j, 1]]
          ],
          Function[j,
            Quiet[Det[{c1Vec, c2Vec,
                       cCenters[[lastChosen]], cCenters[[j]]}]]
          ]
        ];

        (* 選択規則 *)
        new = Which[

          rule === "RZ",
          If[Length[cans] >= 1,
            (* RZ 主: max w，tie -> max geo *)
            wvals = cCenters[[rCans, 4]];
            maxw  = Max[wvals];
            tied  = Select[rCans, cCenters[[#, 4]] >= maxw - eps &];
            If[Length[tied] == 1, tied[[1]],
              gvals = geoScore[#] & /@ tied;
              tied[[ First[Ordering[-gvals]] ]]
            ],
            (* RZ fallback: min w，tie -> max |geo| *)
            wvals = cCenters[[rCans, 4]];
            minw  = Min[wvals];
            tied  = Select[rCans, cCenters[[#, 4]] <= minw + eps &];
            If[Length[tied] == 1, tied[[1]],
              gvals = Abs[geoScore[#]] & /@ tied;
              tied[[ First[Ordering[-gvals]] ]]
            ]
          ],

          rule === "RS",
          If[Length[cans] >= 1,
            (* RS 主: max geo，tie -> max w *)
            gvals = geoScore[#] & /@ rCans;
            maxg  = Max[gvals];
            tied  = Pick[rCans, Map[# >= maxg - eps &, gvals]];
            If[Length[tied] == 1, tied[[1]],
              tied[[ First[Ordering[-cCenters[[tied, 4]]]] ]]
            ],
            (* RS fallback: max |geo|，tie -> min w *)
            gvals   = Abs[geoScore[#]] & /@ rCans;
            maxAbsG = Max[gvals];
            tied    = Pick[rCans, Map[# >= maxAbsG - eps &, gvals]];
            If[Length[tied] == 1, tied[[1]],
              tied[[ First[Ordering[cCenters[[tied, 4]]]] ]]
            ]
          ],

          True,
          rCans[[1]]
        ];

        (* 前進前に状態を保存 *)
        savedNbrs  = neighbors;
        savedLast  = lastChosen;
        hasSaved   = True;

        AppendTo[order, new];
        neighbors  = Complement[#, {new}] & /@ neighbors;
        lastChosen = new;
        justBT     = False
      ]
    ];

    {order, Length[order] == nc}
  ];


(* ── peeling4Df4 と同じインターフェース ─────────────────────── *)
peeling4Df41back[origVers_, faces_, cells_, top_, rule_:"RZ"] :=
  Module[
    {wHat, cc1, localVers, cCenters, neighbors, neighborsInit,
     cansInit, c1Vec, nc, cellVerts, angle4D},

    wHat = {0, 0, 0, 1};
    nc   = Length[cells];

    cellVerts[i_] := Union @@ faces[[ cells[[i]] ]];

    (* 重心を原点へ平行移動 *)
    localVers = (# - Mean[origVers]) & /@ origVers;

    (* C1 を +w 方向に整列（Face-up） *)
    cc1 = Mean[ localVers[[ cellVerts[top] ]] ];
    localVers = Which[
      Abs[(cc1 / Norm[cc1]) . wHat - 1] < 0.0001,
        localVers,
      Abs[(cc1 / Norm[cc1]) . wHat + 1] < 0.0001,
        -localVers,
      True,
        angle4D = ArcCos[Clip[cc1 . wHat / Norm[cc1], {-1, 1}]];
        RotationMatrix[angle4D, {cc1 / Norm[cc1], wHat}] . # & /@ localVers
    ];

    (* セル重心（Face-up 後） *)
    cCenters = Mean[ localVers[[ cellVerts[#] ]] ] & /@ Range[nc];
    c1Vec    = cCenters[[top]];

    (* 隣接リスト（top を除去） *)
    neighbors = Table[
      Select[Range[nc], Function[j,
        j != i && Length[Intersection[cells[[i]], cells[[j]]]] >= 1]],
      {i, nc}
    ];
    neighbors = neighborsInit = Complement[#, {top}] & /@ neighbors;

    (* C2 候補 = C1 の隣接セル *)
    cansInit = neighbors[[top]];

    (* C2 ごとに 1段バックトラックを実行 *)
    Table[
      Module[{c2Vec, res},
        c2Vec = cCenters[[ cansInit[[i]] ]];
        res   = p4f41backPair[
                  c1Vec, c2Vec, cCenters, nc,
                  neighborsInit, top, cansInit[[i]], rule
                ];
        {localVers, res[[1]], res[[2]]}
      ],
      {i, Length[cansInit]}
    ]
  ];

Print["peeling4Df41back loaded."];
Print["  Left condition: Det[{c1,c2,c_k,c_j}] >= -eps  (global, 4D volume form)"];
Print["  k=2: all pass; geo-score = c2_x*cj_y - c2_y*cj_x"];
Print["  RZ: max-w / max-geo  |  RS: max-geo / max-w  |  fallback inverted"];
Print["  Backtrack: 1-step, no consecutive (justBT flag)"];
