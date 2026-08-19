(* ================================================================
   check_5cell_congruence.m
   ================================================================
   正5胞体の20個のSTLが幾何的に同一（合同）かどうかを検証する．

   方針:
     各STLから頂点座標を取り出し，「重心を原点にそろえた上での
     全頂点ペア距離のソート済みリスト」を計算する．
     これは剛体変換（回転・平行移動）不変な等長写像不変量であり，
     20個すべてが一致すれば合同と判断できる．
   ================================================================ *)

stlDir = NotebookDirectory[] <> "STLs/";

(* --- ステップ1: 全STLを読み込み，ユニーク頂点を取得 --- *)
Print["=== 正5胞体 3D展開の合同性検証 ==="];
Print["STLフォルダ: ", stlDir];

files = Sort @ FileNames["unfold_5cell_*.stl", stlDir];
Print["ファイル数: ", Length[files]];

(* 各STLから頂点座標リストを抽出する関数 *)
getUniqueVerts[file_] := Module[{gr, coords},
  gr = Import[file, "Graphics3D"];
  (* Graphics3D から全三角形の頂点を収集 *)
  coords = Cases[gr, Polygon[pts_] :> pts, Infinity] // Flatten[#, 1] &;
  (* 数値丸めてユニーク化（STLの精度限界対策）*)
  DeleteDuplicates[Round[N[coords], 10^-5]]
];

(* --- ステップ2: 各ファイルの不変量（距離リスト）を計算 --- *)
Print["\n--- 各STLの頂点数と距離リスト（先頭10個）---"];

invariants = Table[
  verts = getUniqueVerts[files[[i]]];
  (* 全頂点ペアの距離をソートしたリスト *)
  dists = Sort[
    Table[Norm[verts[[j]] - verts[[k]]],
          {j, Length[verts]}, {k, j+1, Length[verts]}] // Flatten
  ];
  distsRounded = Round[dists, 10^-4];
  If[i <= 5 || i == Length[files],
    Print["  [", i, "] 頂点数=", Length[verts],
          "  距離リスト先頭10: ", Take[distsRounded, 10]]
  ];
  distsRounded,
  {i, Length[files]}
];

(* --- ステップ3: 全ペアが同一かチェック --- *)
Print["\n--- 同一性判定 ---"];
ref = invariants[[1]];
allSame = AllTrue[invariants, # === ref &];

Print["全20個が同一の距離リストを持つか: ", allSame];

If[!allSame,
  (* どれが違うか特定 *)
  diffIdx = Select[Range[Length[files]], invariants[[#]] =!= ref &];
  Print["  異なるSTL: ", diffIdx];
  Print["  差異の例:"];
  Do[
    Print["    [1] vs [", idx, "] 距離リスト長: ",
          Length[ref], " vs ", Length[invariants[[idx]]]];
    diffs = Select[Range[Min[Length[ref], Length[invariants[[idx]]]]],
                   Abs[ref[[#]] - invariants[[idx, #]]] > 0.001 &];
    If[Length[diffs] > 0,
      Print["    最初の不一致 index=", diffs[[1]],
            "  ref=", ref[[diffs[[1]]]],
            "  other=", invariants[[idx, diffs[[1]]]]]
    ],
    {idx, diffIdx}
  ]
];

(* --- ステップ4: 頂点数・辺長の統計まとめ --- *)
Print["\n--- 形状統計（全STL共通であるはず）---"];
verts1 = getUniqueVerts[files[[1]]];
dists1 = Sort[Table[Norm[verts1[[j]] - verts1[[k]]],
                    {j, Length[verts1]}, {k, j+1, Length[verts1]}] // Flatten];
Print["  頂点数（重複除去後）: ", Length[verts1]];
Print["  頂点ペア数: ", Length[dists1]];
Print["  辺長のユニーク値（0.001丸め）: ",
      DeleteDuplicates[Round[dists1, 10^-3]]];
Print["  最小距離: ", N[Min[dists1], 6]];
Print["  最大距離: ", N[Max[dists1], 6]];

(* --- ステップ5: 重心と慣性テンソルで確認 --- *)
Print["\n--- 重心・慣性テンソル主値（各STL）---"];
Table[
  v = getUniqueVerts[files[[i]]];
  ctr = Mean[v];
  vc = v - ConstantArray[ctr, Length[v]];
  inertia = Transpose[vc] . vc;
  ev = Sort[N[Eigenvalues[inertia], 6]];
  If[i <= 3 || i == Length[files],
    Print["  [", i, "] 固有値: ", Round[ev, 10^-3]]
  ],
  {i, Length[files]}
];

Print["\n=== 検証完了 ==="];
