(* ================================================================
   peeling1back.m
   ================================================================
   1段バックトラック付き apple-peel peeling

   アルゴリズム:
     1. グリーディに次の胞を選ぶ（左回り優先, w座標優先）
     2. 候補がなくなって詰まったとき:
          → 直前の1胞だけ取り消して別の候補を試す
     3. 取り消し直後にまた詰まった場合は終了
        （2段以上は戻らない）
     4. 前進できたら justBacktracked フラグをリセット
        → 再びバックトラック可能

   使い方:
     SetDirectory[NotebookDirectory[]];
     Get["peeling1back.m"];
     Get["unfold3DExport.m"];    (* 3D展開・STL出力 *)

     (* 全多胞体を処理 *)
     results = p1bProcessAll[];

     (* 個別処理 *)
     res5 = p1bProcessOne[5];
   ================================================================ *)


(* ================================================================
   Section 1: 補助関数
   ================================================================ *)

(* 胞の頂点番号リスト (1-indexed) *)
p1bCellVerts[ci_, cells_, faces_] :=
  Sort[Union @@ faces[[ cells[[ci]] ]]];

(* 胞隣接リストを構築 *)
p1bBuildAdj[cells_] :=
  Module[{nc, adj},
    nc  = Length[cells];
    adj = Table[
      Sort[Select[Range[nc],
        Function[j, j != i && Length[Intersection[cells[[i]], cells[[j]]]] > 0]]],
      {i, nc}
    ];
    adj
  ];

(* 胞重心を計算 (top をw軸に整列済みの頂点座標 lv を使う) *)
p1bCellCentroids[lv_, cells_, faces_, nc_] :=
  Table[Mean[lv[[ p1bCellVerts[i, cells, faces] ]]],
        {i, nc}];

(* top をw軸（+w方向）に整列させる回転を lv に適用 *)
p1bAlignTopToW[lv0_, cells_, faces_, top_] :=
  Module[{lv, cc1, cc1n, u, wVec, dot, vVec, a, R4},
    lv   = N[lv0 - ConstantArray[Mean[lv0], Length[lv0]]];
    cc1  = Mean[lv[[ p1bCellVerts[top, cells, faces] ]]];
    cc1n = Norm[cc1];
    wVec = {0., 0., 0., 1.};

    If[cc1n > 10^-10,
      u   = cc1 / cc1n;
      dot = Clip[u . wVec, {-1., 1.}];
      Which[
        Abs[dot - 1] <= 10^-4,
          lv,   (* すでに整列済み *)
        Abs[dot + 1] <= 10^-4,
          -lv,
        True,
          vVec = Normalize[wVec - (wVec . u) u];
          a    = ArcCos[dot];
          R4   = (IdentityMatrix[4]
                  + Sin[a]  (Outer[Times, vVec, u] - Outer[Times, u, vVec])
                  + (Cos[a] - 1) (Outer[Times, u, u] + Outer[Times, vVec, vVec]));
          lv . Transpose[R4]
      ],
      lv
    ]
  ];


(* ================================================================
   Section 2: 1段バックトラック peeling（単一 top・f2 版）— v2 修正版
   ================================================================ *)

(* cc    : 胞重心配列 (nc × 4), w軸整列済み
   adj   : 隣接リスト（1-indexed リストのリスト）
   nc    : 胞数
   top   : 開始胞（1-indexed）
   f2    : 2番目の胞（1-indexed）
   戻り値: {order, completeBool}

   バックトラック規則:
     詰まったら直前の1胞を取り消して tried に追加。
     justBacktracked=True のときにまた詰まれば終了。
     前進できたら justBacktracked=False にリセット。

   v2: nbrs をグリーディ版と完全に一致させ，
       バックトラック時に nbrs を save/restore する。
*)
p1bPeelSingleTopF2[cc_, adj_, nc_, top_, f2_] :=
  Module[{order, inOrder, tried, justBT, last, pos,
          nbrs, nbrsInit, savedNbrs, savedLast, hasSaved,
          triedHere, pendingRaw, candidates, bad, prevPos,
          cLast, left, new},

    order     = {top, f2};
    inOrder   = <| top -> True, f2 -> True |>;
    tried     = <||>;                    (* pos -> {tried cells} *)
    justBT    = False;

    (* nbrs 初期化: グリーディ版と同じ構造 *)
    nbrsInit  = Table[DeleteCases[adj[[i]], top], {i, nc}];
    nbrs      = Table[DeleteCases[nbrsInit[[i]], f2], {i, nc}];

    last       = f2;
    hasSaved   = False;
    savedNbrs  = {};
    savedLast  = 0;

    While[Length[order] < nc,
      pos = Length[order];

      (* グリーディ mutation (1): nbrs[last] を未訪問のみにフィルタ *)
      pendingRaw = Select[nbrs[[last]],
                     Function[j, !KeyExistsQ[inOrder, j]]];
      nbrs[[last]] = pendingRaw;

      (* tried フィルタ *)
      triedHere  = Lookup[tried, pos, {}];
      candidates = Select[pendingRaw,
                     Function[j, !MemberQ[triedHere, j]]];

      If[Length[candidates] == 0,
        (* ── 詰まり ── *)
        If[justBT || Length[order] <= 2,
          Break[];
        ];
        bad     = Last[order];
        order   = Most[order];
        inOrder = KeyDrop[inOrder, bad];

        (* nbrs を直前の前進ステップの状態に復元 *)
        If[hasSaved,
          nbrs     = savedNbrs;
          last     = savedLast;
          hasSaved = False
        ];

        prevPos = Length[order];
        tried[prevPos] = Append[Lookup[tried, prevPos, {}], bad];
        justBT  = True,

        (* ── グリーディ選択 ── *)
        cLast   = cc[[last]];
        left    = Select[candidates,
                    Function[j, cLast[[2]] cc[[j,1]] - cLast[[1]] cc[[j,2]] >= 0]];
        new     = If[Length[left] > 0,
                    left[[ First[Ordering[cc[[left, 4]], -1]] ]],
                    candidates[[ First[Ordering[cc[[candidates, 4]], 1]] ]]
                  ];

        (* mutation (2) の前に nbrs を保存 *)
        savedNbrs = Map[List @@ # &, nbrs];  (* deep copy *)
        savedLast = last;
        hasSaved  = True;

        AppendTo[order, new];
        inOrder[new] = True;

        (* グリーディ mutation (2): new を全 nbrs から除去 *)
        nbrs = Table[DeleteCases[nbrs[[i]], new], {i, nc}];
        last = new;
        justBT = False
      ]
    ];

    {order, Length[order] == nc}
  ];


(* ================================================================
   Section 3: 全 top × F2 を列挙する版
   ================================================================ *)

(* 全 (top, f2) 組み合わせについて p1bPeelSingleTopF2 を実行する.
   戻り値: {{top, f2idx, order, completeBool}, ...}
*)
p1bPeelAll[vers_, faces_, cells_] :=
  Module[{nc, adj, lv0c, results, lv, cc, f2Cands, res, order, complete},
    nc     = Length[cells];
    adj    = p1bBuildAdj[cells];
    lv0c   = N[vers - ConstantArray[Mean[vers], Length[vers]]];
    results = {};

    Do[
      (* top をw軸に整列 *)
      lv = p1bAlignTopToW[lv0c, cells, faces, top];
      cc = p1bCellCentroids[lv, cells, faces, nc];

      (* F2 候補 = top の隣接胞 *)
      f2Cands = adj[[top]];

      Do[
        {order, complete} = p1bPeelSingleTopF2[cc, adj, nc, top, f2Cands[[f2i]]];
        AppendTo[results, {top, f2i, order, complete}],
        {f2i, Length[f2Cands]}
      ],
      {top, nc}
    ];

    results
  ];


(* ================================================================
   Section 4: 1つの正多胞体を処理
   ================================================================ *)

(* n 胞体のデータ読み込み・1段バックトラックpeeling・重複チェック・STL出力
   依存: unfold3DExport.m が読み込まれていること
*)
p1bProcessOne[n_] :=
  Module[{dataDir, fname, raw, vers, faces, cells,
          peelingResults, complete, uniqueOrders, seen,
          nOK, nOV, details, top, f2i, order,
          cellEmbs, ovs, bestOrder, bestEmbs, outfile},

    dataDir = NotebookDirectory[] <> "_4DData/";
    fname   = dataDir <> "f" <> ToString[n] <> ".m";

    If[!FileExistsQ[fname],
      Print["File not found: ", fname];
      Return[$Failed]
    ];

    raw   = Get[fname];
    vers  = N[raw[[1]]];
    faces = raw[[4]];
    cells = raw[[6]];

    Print["\n", StringRepeat["=", 50]];
    Print["  ", n, "-cell  (1-step backtrack)"];
    Print[StringRepeat["=", 50]];
    Print["  Cells = ", Length[cells]];

    (* 1段バックトラック peeling *)
    peelingResults = p1bPeelAll[vers, faces, cells];

    complete = Select[peelingResults, #[[4]] &];
    Print["  Complete orders: ", Length[complete]];

    (* 重複除去 *)
    seen         = <||>;
    uniqueOrders = {};
    Do[
      If[!KeyExistsQ[seen, r[[3]]],
        seen[r[[3]]] = True;
        AppendTo[uniqueOrders, r]
      ],
      {r, complete}
    ];
    Print["  Unique complete orders: ", Length[uniqueOrders]];

    If[Length[uniqueOrders] == 0,
      Print["  -> No complete order found."];
      Return[<|"n" -> n, "nOK" -> 0, "nOverlap" -> 0, "details" -> {}|>]
    ];

    (* 重複チェック *)
    nOK      = 0;
    nOV      = 0;
    details  = {};
    bestOrder = Null;
    bestEmbs  = Null;

    Do[
      top   = uniqueOrders[[i, 1]];
      f2i   = uniqueOrders[[i, 2]];
      order = uniqueOrders[[i, 3]];

      cellEmbs = unfoldTo3D[vers, faces, cells, order];
      ovs      = checkOverlaps[order, cells, faces, cellEmbs];

      If[Length[ovs] == 0,
        nOK++;
        AppendTo[details, <|"top" -> top, "f2" -> f2i,
                             "order" -> order, "overlap" -> False|>];
        If[bestOrder === Null,
          bestOrder = order;
          bestEmbs  = cellEmbs
        ];
        Print["  [", i, "/", Length[uniqueOrders], "] top=", top,
              " f2=", f2i, ": OK"],

        nOV++;
        AppendTo[details, <|"top" -> top, "f2" -> f2i,
                             "order" -> order, "overlap" -> True,
                             "nPairs" -> Length[ovs]|>]
      ],
      {i, Length[uniqueOrders]}
    ];

    Print["\n  Summary: ", nOK, " OK, ", nOV, " OVERLAP"];

    (* STL 出力（重複なし最初の結果） *)
    If[bestOrder =!= Null,
      outfile = NotebookDirectory[] <>
                "unfold_" <> ToString[n] <> "cell_1back.stl";
      exportUnfoldSTL[outfile, bestOrder, vers, faces, cells, bestEmbs];
      Print["  -> ", outfile]
    ];

    <|"n"        -> n,
      "nOK"      -> nOK,
      "nOverlap" -> nOV,
      "details"  -> details|>
  ];


(* ================================================================
   Section 5: 全多胞体を一括処理
   ================================================================ *)

p1bProcessAll[] :=
  Module[{ns, results},
    ns      = {5, 8, 16, 24, 120, 600};
    results = Association[];

    Do[
      results[n] = p1bProcessOne[n],
      {n, ns}
    ];

    Print["\n", StringRepeat["=", 50]];
    Print["  GRAND TOTAL (1-step backtrack)"];
    Print[StringRepeat["=", 50]];
    Do[
      If[KeyExistsQ[results, n] && results[n] =!= $Failed,
        Print["  ", n, "-cell: ",
              results[n]["nOK"], " OK, ",
              results[n]["nOverlap"], " OVERLAP"]
      ],
      {n, ns}
    ];

    results
  ];
