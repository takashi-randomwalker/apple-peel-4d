# Wolfram Language ファイル一覧

4D Apple-Peel 展開プロジェクトで使用する Wolfram 言語（`.m`）ファイルの一覧。

---

## アルゴリズム本体

### `peeling4Df.m`
4次元正多胞体に対する Apple-Peel ペーリングの主実装（現行版）。

- **左条件**: xy 2成分外積 $c_{k,y}\,c'_x - c_{k,x}\,c'_y \ge 0$（ステップ $k$ によらず一定）
- **選択規則**: 規則3（左候補から max $w$，なければ min $w$）
- **引数**: `origVers, faces, cells, top`
- **戻り値**: F2 候補ごとに `{rotatedVers, order, successBool}` のリスト

---

### `peeling4Df3.m`
左条件を xyz 3成分行列式に改良した実験版（未テスト）。

- **左条件**: $k \ge 2$ で $\det(\tilde{c}_{k-1},\,\tilde{c}_k,\,\tilde{c}_j) \ge 0$（xyz 射影），$k = 1$ は xy 2成分にフォールバック（直前セル未定義のため）
- **選択規則**: 規則3（`peeling4Df.m` と同一）
- **引数・戻り値**: `peeling4Df` と同一

---

### `peeling4Df_rules.m`
R1 / R2 / R3 の三規則すべてに対応した 4D ペーリング実装。

- **左条件**: xy 2成分外積（`peeling4Df.m` と同一）
- **選択規則**:
  - `"maxphi"` (R1): 左候補から max $\varphi$（最大左旋回角）
  - `"loxo"` (R2): 左候補から min $\varphi$（loxodrome 近似）
  - `"maxw"` (R3): 左候補から max $w$（剥き軸優先）
- **Darboux フレーム**: xy 射影版（$\hat{\mathbf{f}}^{xy}$, $\hat{\mathbf{l}}^{xy}$）を定義
- **補助関数**: `p4rXYFrame`, `p4rAngleFromForward`

---

### `peeling3DLoxo.m`
3次元多面体に対する Apple-Peel ペーリング（R1 / R2 / R3 対応）。

- **左条件**: $\det(\mathbf{c}_{k-1}, \mathbf{c}_k, \mathbf{c}_j) \ge 0$（3成分行列式，全ステップ厳密版）
- **選択規則**:
  - `"maxphi"` (R1): 左候補から max $\varphi$
  - `"loxo"` (R2): 左候補から min $\varphi$（loxodrome 近似）
  - `"maxz"` (R3): 左候補から max $z$（緯度優先）
- **主要関数**:
  - `p3lDarbouxFrame[cPrev, cCurr]` → `{nHat, fHat, lHat}`
  - `p3lAngleFromForward[fHat, lHat, cj]` → $\varphi_j$
  - `p3lPeelPair[cc, adj, nf, top, f2, rule]` → `{order, completeBool}`
  - `p3lRunAll[vers, faces, rule]` → 全 $(F_1, F_2)$ ペアの結果
  - `p3lFromBuiltin[name, rule]` → `PolyhedronData` から読み込んで実行
  - `p3lCompareRules[]` → 全 Platonic 立体で規則2・3を比較

---

### `peeling1back.m`
1段バックトラック付き 4D ペーリング。

- グリーディ選択で詰まった際に直前の1胞だけ取り消して別の候補を試す
- 取り消し直後にも詰まった場合は終了（2段以上は戻らない）
- 前進できたら `justBacktracked` フラグをリセットし再びバックトラック可能
- **主要関数**:
  - `p1bProcessOne[n]` → $n$-胞体を処理
  - `p1bProcessAll[]` → 全多胞体を処理

---

## 実行スクリプト

### `run4DPeeling.m`
`peeling4Df.m` を全6種の正多胞体に適用するスクリプト。

- 対象: 5 / 8 / 16 / 24 / 120 / 600-cell
- 結果を `Perfect / Possible / Impossible` に分類して出力
- 結果を変数 `results4D` に格納（`ans4D.mx` として保存するオプション付き）

---

### `run4DPeelingR3.m`
`peeling4Df3.m`（xyz 行列式左条件）を全6種の正多胞体に適用するスクリプト。

- 結果を `results4DR3` に格納（`ans4DR3.mx` として保存するオプション付き）
- **未実行**

---

### `run4DPeelingAllRules.m`
`peeling4Df_rules.m` を使い R1 / R2 / R3 の三規則を全6種の正多胞体に一括実行するスクリプト。

- 成功数・ユニーク成功数を規則ごとに比較表として出力
- 結果を `rulesResults4D` に格納（`rules_results_4D.mx` として保存するオプション付き）

---

### `compare_TI_R1_R3.m`
切頂二十面体（Truncated Icosahedron）における R1 と R3 の詳細比較スクリプト。

- `peeling3DLoxo.m` を流用し R1 (`maxphi`) ラッパーを追加定義
- 面種別 (pentagon / hexagon) の選択順列を出力し，zigzag（R1）vs spiral（R3）のパターンを比較

---

## 展開・検証・出力

### `unfold3DExport.m`
成功した 4D ペーリング結果を 3D 空間に実現し STL 出力・重なり判定を行う。

- **SVD 射影**: 各セルの4D頂点を最良3D超平面に等長射影（`unfoldCellLocal3D`）
- **Procrustes 整列**: 共有面頂点の RMSD を最小化する回転・平行移動を算出（`unfoldProcrustes`）；鏡映補正付き
- **SAT 重なり判定**: 面法線・辺外積を候補軸とする Separating Axis Theorem（`satOverlap`, `checkOverlaps`）；AABB 事前フィルタ付き
- **STL 出力**: `exportUnfoldSTL[vers, faces, cells, order, filename]`
- **可視化**: `showUnfolding[vers, faces, cells, order]` → `Graphics3D`
- **一括処理**:
  - `processOnePolytope[n]` → $n$-胞体の全成功順序を処理
  - `processAllResults[]` → 全多胞体を処理

---

## 解析・確認スクリプト

### `check_fallback.m`
fallback（左候補なし時の min-$w$ 選択）の発火状況と，fallback なし版との成功数比較を行うスクリプト。

- **`peeling4DfLog`**: fallback 発火カウント付き版
- **`peeling4DfNoFB`**: fallback なし版（左候補なし → 即 Break）
- 5 / 8 / 16 / 24-cell を対象に両版の成功数・fallback 回数を比較
- 3D 正十二面体の R3 成功数も確認

---

### `check_overlap.m`
`ans4D.mx` の全成功展開図について SAT 重なり判定を実行し，valid net（重なりなし）数を集計するスクリプト。

- `peeling4Df.m` と `unfold3DExport.m` を読み込んで使用
- 5 / 8 / 16 / 24 / 120-cell を対象
- valid / overlap 数と割合を出力

---

## データファイル（`_4DData/`）

各ファイルは `{vers, edges, vertexNeighbors, faces, faceNeighbors, cells}` の6要素リストを返す。
`vers` は4D頂点座標，`faces` は頂点インデックスリスト，`cells` は面インデックスリスト。

| ファイル | 多胞体 | セル数 | ファイルサイズ |
|----------|--------|-------:|---------------:|
| `f5.m`   | 5-cell（正五胞体） | 5 | 0.7 KB |
| `f8.m`   | 8-cell（正八胞体・超立方体） | 8 | 2.5 KB |
| `f16.m`  | 16-cell（正十六胞体） | 16 | 2.3 KB |
| `f24.m`  | 24-cell（正二十四胞体） | 24 | 6.6 KB |
| `f120.m` | 120-cell（正百二十胞体） | 120 | 119 KB |
| `f600.m` | 600-cell（正六百胞体） | 600 | 127 KB |

---

## 計算済みデータ

| ファイル | 内容 |
|----------|------|
| `ans4D.mx` | `results4D`：`peeling4Df.m` による全6種の計算結果（Mathematica バイナリ） |
| `archimedean_results.mx` | Archimedean 13種の R1 / R2 / R3 比較結果 |
| `221110ArchimedeanRightHand.m` | Archimedean 立体の右手座標系ペーリング結果（数値データ） |

---

## 依存関係

```
peeling4Df.m          ←─ run4DPeeling.m
peeling4Df3.m         ←─ run4DPeelingR3.m
peeling4Df_rules.m    ←─ run4DPeelingAllRules.m
peeling3DLoxo.m       ←─ compare_TI_R1_R3.m
peeling4Df.m          ←─ unfold3DExport.m  ←─ peeling1back.m
peeling4Df.m          ←─ check_fallback.m
peeling4Df.m          ←─ check_overlap.m
unfold3DExport.m      ←─ check_overlap.m
_4DData/fN.m          ←─ （上記すべての実行スクリプト）
```
