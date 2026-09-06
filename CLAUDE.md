# CLAUDE.md — 4D Apple-Peel 展開プロジェクト

## プロジェクト概要

3次元多面体の **Apple-Peel 展開**（Yoshino, Chaidee）を **4次元正多胞体**に拡張するプロジェクト。

- 実装言語：**Mathematica**（`.m` スクリプト + `.nb` Notebook）
- ドキュメント：LaTeX / Markdown
- 対象：6種の正多胞体（5胞体・8胞体・16胞体・24胞体・120胞体・600胞体）

詳細な計算ログ・廃止された実装の比較・セッション別修正履歴は **`HISTORY.md`** を参照。
アルゴリズムの理論的詳細・等変性証明は **`summary.tex` / `summary_en.tex`** を参照。

---

## ディレクトリ構成（2026-08-19 に整理）

トップレベルは342ファイルまで膨れていたので話題ごとに分けた。**現在トップレベルにあるのは論文（`.tex`/`.md`）と図（`.pdf`/`.png`）だけ**で90ファイル。

```
src/       スクリプト85本を8つに分類
  peel3d/    13  3D apple-peel（本体・実行・可視化）
  peel4d/    30  4D apple-peel（本体・実行・展開・分析）
  bfsnet/     8  face-rotation BFS ネット（bfs_sat.m を含む）
  uniform/   14  Wythoff 生成・一様多胞体の実行とエクスポート
  symmetry/   9  対称性の上限・扇形ネット・曲率（2026-08-19 の仕事）
  tiebreak/   6  タイブレーク規則の調査
  random/     3  ランダム多面体
  figures/    2  論文・発表用の図の生成
data/      *.mx 34件（計算結果。.gitignore 対象なので git には入らない）
notebooks/ *.nb 13件
attic/     oneoff/ 35（使い捨ての check_*/test_*/gen_*/debug_*）
           outputs/ 19（loose STL・progress ログ）
           figsrc/ 12（図の Illustrator/EPS 原本）
           misc/ 18（.docx・.mov・スクショ・投稿zip 等）
_4DData/  _4DData_uniform/         多胞体の組合せデータ（**移動していない**）
uniform_nets/  STLs_v4/  face_rotation_nets/   公開用ネット（**移動していない**）
```

**移動しなかったものの理由**：`_4DData/` 系は参照が多く公開データでもある。`uniform_nets/` 等は **`.gitignore` に `!uniform_nets/stl/*.stl` のようなパス依存の例外規則があり，動かすと STL が除外されてしまう**。論文と図をトップレベルに残したのは，JoCG 投稿直前に `\includegraphics` のパスを触るリスクを避けるため（投稿後に `papers/` へ移す予定）。

**スクリプト間の参照は `baseDir <> "src/<topic>/<name>"` の形に書き換え済み**（50箇所）。`.mx` の参照も `data/` 付きに書き換え済み（46箇所）。`baseDir` はプロジェクトルートのままなので，`baseDir` の定義自体は変更していない。

> 移動後の検証：全42箇所の参照先が実在することを確認し，`src/symmetry/symmetry_ceiling.wls` を新しい場所から実走して移動前と同一の結果（24/8/1/6）を確認。`face_rotation_net_paper.tex` と `jcdcg3_slides.tex` も 12ページでビルド通過。

---

## アルゴリズム（現行版）

### 3次元版（`src/peel3d/peeling3DLoxo.m`）

#### 前処理
- **Face-up 回転**：`p3lAlignTopToZ[vers, faces, top]` で開始面 F1 の重心を +z 軸に揃える
- **全ペア評価**：`p3lRunAll[vers, faces, rule]` が全 (F1, F2) ペアを評価

#### 1ステップの選択
1. `last` の未訪問隣接面集合を `pending` とする
2. `pending` が空なら終了（失敗）
3. **大域的右半空間条件**で候補を絞る（単一候補でも適用）：
   - `Det[{c_1, c_k, c_j}] <= eps`（変数名は `leftCands`）
   - `c_1 = cc[[top]]`（開始面の重心，固定参照）；`eps = 10^-10`（数値誤差吸収）
   - 右半空間 = 外側から見て CW 方向の候補
4. `rightCands` が空でなければ規則 r で選択；空ならフォールバック（pending 全体から）

**大域参照の理由**：c_{k-1}（局所参照）は各ステップで変化し螺旋方向の一貫性を保証できない。c_1（固定参照）で全ステップの大域的螺旋方向を保証。

#### 選択規則（Det ベース統一，2026-05-12）

| 規則 | 右候補あり | フォールバック |
|------|-----------|--------------|
| **RS** = R1（Spiral） | min `Det[{c_1,c_k,c_j}]`（最大右旋回） | Darboux frame の最小 φ |
| **RZ** = R3（Zonal） | max z；同値（差 ≦ 10⁻¹⁰）なら min Det | min z；同値なら min Det |
| R2（loxodrome） | max Det | max \|Det\| |

> **⚠️ 論文方針（2026-05-10）**：`paper_draft.tex` では **R2 は掲載しない**（記録としてのみ CLAUDE.md・summary.tex に保持）。

**RZ 二次基準の理由**：正十二面体のように高対称な多面体では Face-up 後に複数候補が理論上同 z。浮動小数点誤差で ~10⁻¹⁶ の差が生じリスト順タイブレークすると等変性が崩れる。min Det で等変性回復（旧版 min φ と実質等価）。

Darboux フレーム：`p3lDarbouxFrame[cPrev, cCurr]` → `{nHat, fHat, lHat}`，`p3lAngleFromForward` で φ 計算（RS フォールバックのみ使用）。

#### フォールバックなし版
`src/peel3d/run_platonic_updated.m` の `p3lPeelPairNoFB2`：ステップ 4 でフォールバック使わず即終了（右条件は適用）。

### 4次元版（`src/peel4d/peeling4Df4.m`，現行・推奨）

- **剥き軸**：w 軸（Face-up 配向）
- **左条件**：`Det[{c_1, c_2, c_k, c_j}] >= -eps`（4次元体積形式，**c1–c2 平面グローバル参照**）
  - 3D 版 `Det[{c_1, c_k, c_j}]` の直接類推
  - Face-up 後 c_1 = (0,0,0,w_1) → `-w_1 * Det_xyz[{c_2, c_k, c_j}] >= -eps` と等価
  - k=2 のとき Det[{c1,c2,c2,c_j}] = 0 → 全候補通過（特別扱い不要）
  - **等変性**：A ∈ SO(4) は det(A)=1 → `Det[{Ac1,Ac2,Ac_k,Ac_j}] = Det[{c1,c2,c_k,c_j}]`

- **選択規則 RZ**：左候補から max w；同点（差 ≦ 10⁻¹⁰）なら geo-score タイブレーク；フォールバック min w，同点なら max |geo-score|
- **選択規則 RS**：左候補から max geo-score；同点なら max w；フォールバック max |geo-score|，同点なら min w

#### geo-score（タイブレーカー）

| k | 式 | 備考 |
|---|-----|------|
| k=2 | `c2_x·cj_y − c2_y·cj_x`（xy クロス積） | +w 軸から見た c2→cj 反時計回り角度。k=2 は c_k=c_2 で Det≡0 となるため |
| k≥3 | `Det[{c1,c2,c_k,c_j}]` | 全候補が ≈0（Det=0 縮退）なら xy クロス積にフォールバック（2026-05-12 追加） |

**k≥3 Det=0 フォールバック**：8胞体（セル重心がすべて ±e_i）と 120胞体南極付近で頻発するリスト順依存を解消し，120胞体を Perfect に押し上げた重要修正。

> **旧版 `src/peel4d/peeling4Df.m` / `src/peel4d/peeling4Df3.m`** は左条件が異なる（xy 2成分 / xyz ローカル参照）。現役は `src/peel4d/peeling4Df4.m` のみ。論文からも旧版の議論は削除済み。

#### 1段バックトラック（`src/peel4d/peeling4Df4_1back.m`）

詰まったとき1ステップ戻って別候補を試す拡張（justBT フラグで連続バックトラック禁止）。全6胞体で **グリーディ版と完全に同一結果**（改善なし，2026-05-11）。

---

## 主要ファイル（現役）

### 4D 本体・実行スクリプト

| ファイル | 役割 |
|----------|------|
| `src/peel4d/peeling4Df4.m` | 4D本体（c1–c2 平面グローバル参照 + RZ/RS） |
| `src/peel4d/peeling4Df4_1back.m` | 1段バックトラック版（グリーディと同値） |
| `src/peel4d/run4DGlobal_update.m` | 5/8/16/24胞体の v4 計算（k≥3 Det=0 fix 適用） |
| `src/peel4d/run4DGlobal_120cell.m` | 120胞体の v4 計算 |
| `src/peel4d/run4DGlobal_nofallback.m` | RZ/RS × with/without fallback 比較 |
| `src/peel4d/run4DGlobal_RSRZ.m` | RS vs RZ 比較 |
| `src/peel4d/run4DBacktrack.m` | 1段バックトラック計算 |
| `src/peel4d/run4DExportSTLv4.m` | v4 から valid net を STLs_v4/ に出力 |

### 3D 本体・実行スクリプト

| ファイル | 役割 |
|----------|------|
| `src/peel3d/peeling3DLoxo.m` | 3D本体（RS/RZ + R2） |
| `src/peel3d/run_platonic_updated.m` | 正多面体 5種計算（fallback 有無） |
| `src/peel3d/run_platonic_strict.m` | 厳密条件（eps<0）計算 |
| `src/peel3d/run_archimedean_faceup.m` | アルキメデス 13種（fallback あり） |
| `src/peel3d/run_archimedean_nofallback.m` | アルキメデス 13種（fallback なし） |

### 可視化・展開・STL 出力

| ファイル | 役割 |
|----------|------|
| `src/peel4d/unfold3DExport.m` | SVD射影・Procrustes整列・SAT判定・STL出力 |
| `src/peel4d/unfold4D.m` | 4D展開図の3D可視化（`data/ans4DGlobal_v4.mx` + `src/peel4d/unfold3DExport.m`） |
| `src/peel3d/visualize_platonic_nets_v2.m` | 正多面体展開図描画 |
| `src/peel3d/visualize_archimedean_nets.m` | アルキメデス多面体展開図描画 |
| `gen_120cell_faceup_nets.m` | 120胞体 cell-up vs face-up 比較（v4 fix 適用） |

### データファイル（現役）

| ファイル | 内容 |
|----------|------|
| `data/ans4DGlobal_v4.mx` | **現行最新**：4D 全6胞体の計算結果（v4） |
| `data/ans4DGlobal_nofb.mx` | フォールバックなし比較（サマリのみ） |
| `data/ans4DGlobal1back.mx` | 1段バックトラック結果 |
| `data/dataPlatonic.mx` | 正多面体 × R1/R2/R3 × fallback 有無 |
| `data/dataPlatonic_strict.mx` | 正多面体 × RS/RZ × fallback 有無（eps<0） |
| `data/archimedean_faceup_results.mx` | アルキメデス（fallback あり，Det ベース） |
| `data/archimedean_nofallback_results.mx` | アルキメデス（fallback なし，Det ベース） |
| `_4DData/f5.m` 〜 `f600.m` | 正多胞体データ（頂点・面・セル） |
| `STLs_v4/` | 280枚の valid net STL（重なりなし保証） |

### ドキュメント

| ファイル | 内容 |
|----------|------|
| `paper_draft.tex` / `.pdf` | 投稿用論文ドラフト |
| `summary.tex` / `.pdf` | アルゴリズム詳細（最も包括的，日本語） |
| `summary_en.tex` / `.pdf` | 英語版 summary |
| `archimedean_results.tex` / `.pdf` | Archimedean 比較結果 |
| `literature_review_spiral_unfolding.md` / `.tex` / `.pdf` | 既往研究レビュー |
| `commentsOpus260516.md` | Opus 4.7 レビューコメント |
| `HISTORY.md` | 計算ログ・廃止実装・セッション履歴 |

> 旧版（`src/peel4d/peeling4Df.m`, `src/peel4d/peeling4Df3.m`, `peeling1back.m`, `data/ans4DGlobal.mx`, `data/ans4DGlobal_v2.mx`, `data/ans4DGlobal_v3.mx`, `STLs/`, `STLs_v3/` など）は `HISTORY.md` 参照。

---

## 計算結果（最新確定）

### 4D：`src/peel4d/peeling4Df4.m`（v4，k≥3 Det=0 fix，2026-05-12）

結果：`data/ans4DGlobal_v4.mx`

| 多胞体 | (C1,C2)総数 | True | Unique | 分類 |
|--------|:-----------:|-----:|-------:|------|
| 5胞体  | 20          | 20   | 20     | **Perfect** |
| 8胞体  | 48          | 48   | 48 (geo-unique=1) | **Perfect** |
| 16胞体 | 64          | 20   | 20     | Possible |
| 24胞体 | 192         | 192  | 192    | **Perfect** |
| 120胞体 | 1,440      | 1,440| 1,440  | **Perfect** |
| 600胞体 | 2,400      | 0    | 0      | **Impossible** |

### Valid net 比率（`STLs_v4/`，2026-05-12）

| 多胞体 | Unique | Valid | Valid% |
|--------|:------:|------:|-------:|
| 5/8/16/24胞体 | 280  | 280   | 100% |
| 120胞体 | 1,440  | 0     | 0%     |
| **合計** | **1,720** | **280** | **16.3%** |

- 5/8/16/24胞体：全 order が 3D 印刷可能な valid net
- 120胞体：全 order に自己交差あり（Perfect な ordering でも valid net は生成しない）
- 600胞体：ordering 自体が 0

### フォールバックなし比較（`src/peel4d/peeling4Df4.m` noFB）

| 多胞体 | RZ(w) | RZ(n) | RS(w) | RS(n) |
|--------|------:|------:|------:|------:|
| 5/8/24胞体 | Perfect | Perfect | Perfect | Perfect |
| 16胞体 | Possible/20 | Possible/20 | Impossible | Impossible |
| 120胞体| **Perfect/1440** | Possible/978 | Impossible | Impossible |
| 600胞体| Impossible | Impossible | Impossible | Impossible |

- フォールバックが有効な唯一の組み合わせは **120胞体 × RZ**（462ペア=32%がフォールバック経由）
- 4次元では RZ が RS を大幅に上回る（RS は 16胞体・120胞体で Impossible）

### 3D 正多面体（`data/dataPlatonic.mx`，2026-05-06）

全結果が 0% か 100%（面推移性による等変性理論と一致）。
- with fallback：全 5 種・全規則で 100%
- no fallback：Tetrahedron のみ全規則で 100%；他は R2 のみ 0%（他規則は 100%）
- 厳密条件（eps<0）：Tetrahedron のみ no-fallback で 100%，他は 0%

### 3D アルキメデス（`data/archimedean_faceup_results.mx`，2026-05-12 Det ベース）

| 多面体 | F | RS(w) | RZ(w) | 備考 |
|--------|:-:|------:|------:|------|
| TruncatedTetrahedron | 8 | 66.7% | 33.3% | RS > RZ の唯一例 |
| TruncatedOctahedron | 14 | 41.7% | 100% | RZ > RS |
| TruncatedCuboctahedron | 26 | 0% | 100% | RZ > RS |
| TruncatedIcosahedron | 32 | 0% | 100% | RZ > RS |
| TruncatedIcosidodecahedron | 62 | 0% | 66.7% | RZ > RS |
| **SnubCube** | 38 | 20% | 40% | RS/RZ とも fallback 効果あり（24→0, 48→24） |
| その他7種 | — | 0% | 0% | Cuboctahedron, TruncatedCube, Rhombicuboctahedron, Icosidodecahedron, TruncatedDodecahedron, Rhombicosidodecahedron, SnubDodecahedron |

- フォールバック効果は SnubCube 1 種のみ；他 12 種は fallback 有無で結果不変
- RZ 結果は Det ベース選択導入前後で完全不変（min-φ ≡ min-Det）
- RS 結果は Det ベースで大幅変化（旧版で 100% だった 3 種が 0% に）
- 鏡像対称性：RZ で全 13 種一致；キラルの SnubCube は鏡像間で成功ペアが入れ替わる

---

## 3D ランダム多面体実験結果（2026-06-11）

球面上のランダム母点 N 個から ConvexHullMesh（Delaunay 三角分割）と DualPolyhedron（Voronoi 双対）を生成し，Apple-Peel VorRZ/DelRZ 成功率を比較。

### 主要スクリプト

| ファイル | 内容 |
|----------|------|
| `src/random/run_random_polyhedra.m` | N=30 の Delaunay vs Voronoi 詳細比較（1 インスタンス + 20 試行） |
| `src/random/run_random_scaling.m` | N = {20,30,50,100,200} × 4 条件のスケーリング実験 |
| `src/random/run_repulsive_comparison.m` | Thomson Coulomb 緩和による均一分布の効果実験（N=30） |

### Delaunay vs Voronoi 比較（N=30，20 試行平均）

| 条件 | RS | RZ |
|------|---:|---:|
| Delaunay（全三角形，3-正則，|F|≈2N-4） | ~0% | ~1% |
| Voronoi 双対（混合多角形，平均次数≈5.6） | ~0% | ~30% |

- Delaunay（3-正則）は 4D 16-cell/600-cell 型（Possible/Impossible）に対応
- Voronoi（6-正則）は 4D 8-cell/120-cell 型（Perfect）に対応
- Euler 制約 Σ(6-k)=12：全 Voronoi 双対で 12 個の次数欠陥（pentagon 等）が存在

### N スケーリング（Voronoi RZ）

| N | Vor RZ |
|--:|-------:|
| 20 | ~51% |
| 30 | ~22% |
| 50 | ~9% |
| 100 | ~0.1% |
| 200 | ~0% |

N 増大とともに単調減少。面次数は N→∞ で 6 に近づくが成功率は 0 に向かう。

### Thomson 緩和実験（2026-06-11）

Thomson Coulomb 勾配降下法（lr=0.005）で母点を球面上で均一化し，VorRZ 成功率を測定。

実験1（N=30，steps={0,10,50,200,500,2000}，10 試行/レベル）：

| steps | VorRZ |
|------:|------:|
| 0 | 27.9% |
| 10 | 30.6% |
| 50 | 26.2% |
| 200 | 25.5% |
| 500 | 21.0% |
| 2000 | 27.0% |

実験2（N={20,30,50}，random vs Thomson steps=2000，各 10 試行）：

| N | random | Thomson |
|--:|-------:|--------:|
| 20 | 44.7% | 43.9% |
| 30 | 32.6% | 22.1% |
| 50 | 5.9% | 6.5% |

**結論：Thomson 緩和は成功率を改善しない。** 差はすべて統計誤差（σ≈5–10%）の範囲内。

### 主結論

**成功率を決めるのは面推移性（対称性の存在）であり，局所的な幾何正則性（点配置の均一さ・面の六角形化）ではない。**

- 正多胞体が Perfect なのは面推移群による等変性（0% か 100% の二択）
- ランダム多面体は対称性がないため Thomson 均一化・N によらず ~0–50%（N が大きいほど低い）
- この結果は論文 Future Work「S² ランダム凸包での scaling 解析」の実証的根拠となる

---

## Face-rotation BFS ネット結果（別論文，2026-06-09）

### アルゴリズム概要

BFS spanning tree を根セルから構築し，各子セルを親との共有面周りに 4D 回転させて同じ 3D 超平面に展開する。Apple-peel とは独立な展開法。

1. **Face-up**：根セルの重心を +w 軸に揃える（`faceUpForRootP`）
2. **BFS 展開**：辺 (親, 子) ごとに `unfoldTFP` で 4D アフィン変換を累積
3. **射影**：w 座標を捨て，残りの xyz が 3D ネット

### 計算結果（`src/bfsnet/face_rotation_net_all4D_v2.m`，2026-06-09）

| 多胞体 | セル数 | 根の数 | Valid | 結果 |
|--------|:------:|:------:|------:|------|
| 5-cell   |   5 |   5 |   5 | ALL VALID |
| 8-cell   |   8 |   8 |   8 | ALL VALID |
| 16-cell  |  16 |  16 |  16 | ALL VALID |
| 24-cell  |  24 |  24 |  24 | ALL VALID |
| 120-cell | 120 | 120 | 120 | ALL VALID |
| 600-cell | 600 | 600 | 600 | ALL VALID |
| **合計** | — | **773** | **773** | **全有効** |

**定理**：全6種の正凸 4-多胞体について，任意の根セルから生成した face-rotation BFS ネットは有効（自己交差なし）。

### 重要なバグと修正

**バグ**：`RegionMember`（Mathematica，閉領域）を用いた頂点帰属検査は，隣接しない 2 セルが辺・頂点のみを共有するだけで（境界接触，体積 = 0）"重複あり" と誤判定する。正四面体・正八面体のように尖ったセルで頻発。

**修正**：体積基準 `RegionMeasure[RegionIntersection[h1, h2]] > 10^-8` に変更。境界接触は無視，真の内部重複のみ検出。

- バグの影響：16-cell 0/16（誤）→ 16/16，24-cell 0/24（誤）→ 24/24，600-cell 0/600（誤）→ 600/600
- D&H (2022) が「16-cell は全スパニングツリーで有効」と証明済みであり整合

### Apple-peel との対比（主要な対照例）

| 多胞体 | Apple-peel 有効順序 | Apple-peel 有効ネット | BFS 有効根 |
|--------|:-------------------:|:--------------------:|:----------:|
| 120-cell | 1440/1440 (100%) | **0/1440 (0%)** | **120/120** |
| 600-cell | **0/2400 (0%)** | — | **600/600** |

有効な ordering の存在と有効な 3D ネットの生成は独立した性質。

### 主要スクリプト

| ファイル | 内容 |
|----------|------|
| `src/bfsnet/face_rotation_net_all4D_v2.m` | 全6胞体の BFS ネット検証（体積基準，**現行推奨**） |
| `src/bfsnet/face_rotation_net_small4D.m` | 5/8/16/24-cell 検証（v1，RegionMember 版，歴史的参考） |
| `src/bfsnet/face_rotation_net_600cell.m` | 600-cell 専用検証スクリプト |
| `check_16cell_detail.m` | バグ調査：RegionMember vs 体積基準の比較 |
| `src/bfsnet/face_rotation_net_120cell.m` | 120-cell BFS ネット可視化（net PNG 出力） |
| `face_rotation_net_viewer.nb` | BFS アニメーション付き可視化 Notebook |

### 論文

`face_rotation_net_paper.tex`（6 ページ，arXiv 投稿可）

- タイトル："Face-Rotation BFS Nets of the Six Regular Convex 4-Polytopes"
- 著者：**Takashi Yoshino 単著**（apple-peel 論文 paper_draft.tex とは別著者構成）
- D&H (2022) の 5/8/16-cell 結果を拡張し，24/120/600-cell の open cases を BFS の範囲で解決
- 未解決問題：(1) 計算に依存しない理論的証明，(2) 全スパニングツリーへの拡張
- セクション構成：1. Introduction / 2. Preliminaries（2.1 Background, 2.2 Face-Rotation BFS Algorithm）/ 3. Results（3.1 Main result, 3.2 Contrast with Apple-Peel Unfolding）/ 4. Open Problems

#### ジャーナル投稿可能性（2026-06-25 評価）

**強み**
- D&H (2022) の未解決3ケース（24/120/600-cell）を BFS 族について解決
- 120-cell（全順序が有効なのに有効なネットが0件）と 600-cell（BFS ネット全根有効なのに apple-peel 順序不可能）の対比が genuine に novel：ordering validity と net validity の独立性を実証
- 体積基準の validity test（RegionMember の誤判定を修正）も技術的貢献

**弱み**
- 証明が純粋な計算検証のみ（なぜ成り立つかの理論的説明なし）
- D&H (2022) が全スパニングツリーで証明したのに対し，本論文は BFS 族のみ
- 6ページと短く，スタンドアローンの journal 論文としては内容が薄く見えうる

**投稿先候補と現実的評価**

| ジャーナル | 可能性 | 備考 |
|-----------|:------:|------|
| Journal of Computational Geometry (JoCG) | 中〜高 | 第1候補．オープンアクセス，計算検証結果を受け入れる実績あり |
| Graphs and Combinatorics (Springer) | 中 | 短い組合せ幾何論文も掲載 |
| Discrete Mathematics | 中 | 幅広い離散数学対象，計算的証明も可 |
| CGTA (Elsevier) | 低〜中 | D&H (2022) の掲載誌で最自然だが，理論的貢献の薄さを指摘される恐れ |

**推奨戦略**
1. arXiv に投稿して先行性を確保（即可能）
2. JCDCG^3 採否通知（2026-07-03）を参考に journal 投稿方針を決定
3. JoCG を第1候補として投稿

**注意**：「計算検証のみ」スタイルは査読者の反応が分かれる。「なぜ成り立つかが不明」として major revision / reject になるリスクはある。ただし未解決問題の解決と独立性の発見という2点は publishable な核になる。

#### arXiv 番号の区別（混同注意）

本プロジェクトには 2 つの arXiv 論文が関わる。混同しないこと。

| arXiv 番号 | 論文 | 役割 |
|------------|------|------|
| **2605.30373** | apple-peel 論文（`paper_draft.tex`） | `face_rotation_net_paper.tex` が `\bibitem{YoshinoChaidee2026}` で参照する companion paper |
| **2604.16204** | 旧 3D apple-peel 論文（`peeling3Df` 使用） | `paper_draft.tex` 内の `\bibitem{Yoshino2026arXiv}` で参照される earlier work |

つまり 2604.16204 は `face_rotation_net_paper.tex` からは**直接参照されない**（apple-peel paper を経由した間接参照）。

#### 2026-06-24 セッションの修正

| 修正 | 内容 |
|------|------|
| ORCID 削除 | Acknowledgements 内の `\paragraph{ORCID.}` を削除（著者 footnote に既出のため重複） |
| Figure 1 caption | "rooted at cell~1" を削除（index 情報は不要） |
| 単著対応 | "The authors declare" → "The author declares" |
| ダングリング ref 修正 | `\label{sec:intersect}` 削除 + 該当 `\ref` 参照削除 |
| Open Problems 段落削除 | 「apple-peel 自己交差は distant layer」記述は実測（隣接帯間が主体，`src/peel4d/analyze_120cell_overlap_bands.m` 由来）と矛盾するため段落ごと削除 |
| §2 タイトル | "Theory" → "Preliminaries"（§2.2 が手順記述のため "Theory" は不適切） |
| Figure 参照追加 | Proof 内に `Figure~\ref{fig:allnets}` 参照を追加（孤立図の解消） |
| 変数 $V$ 改名 | Step 2 の visited set $V$ → $\mathrm{Vis}$（vertex set $V$ との衝突解消） |
| §3 タイトル | "Result" → "Results"，§3.1 "Main Result" → "Main result" |
| Algorithm 1 コメント統一 | `\algorithmiccomment` → `\Comment`（algpseudocode 標準への統一） |
| Algorithm 1 clip 廃止 | `\mathrm{clip}(x,-1,1)` → `\max(-1,\min(1,x))`（記号が読者に通じない懸念のため） |
| US 綴り統一 | organised→organized，neighbour(s)→neighbor(s)（4箇所），initialised→initialized，normalised→normalized（2箇所），normalise→normalize，summarises→summarizes，visualisation→visualization（apple-peel paper が US で書かれているため整合） |
| 文章修正 | "v belongs to r" → "v is a vertex of r"；"As a corollary, the result contrasts sharply with..." → "By contrast, the result diverges sharply from..." |
| STL ファイル追加 | Figure 1 caption と Code and Data Availability 節に STL 公開の旨を追記 |
| GitHub push | `face_rotation_nets/` ディレクトリに全6胞体の STL（root=1）を追加；`.gitignore` に例外規則を追記 |

#### 2026-06-25 セッションの修正（詳細レビュー対応）

| 修正 | 内容 |
|------|------|
| H1: G(P) 定義の誤記修正 | "edges = shared codimension-1 faces" → "edges = shared ridges, i.e., shared 2-faces"（4-多胞体の codimension-1 face は cell 自身であるため誤り） |
| H1': cell degree 記述の同様の誤記修正 | "cell degree denotes... (sharing a codimension-1 face)" → "(sharing a ridge, i.e., a 2-face)"（同上の理由） |
| H2: "non-adjacent" の曖昧さ解消 | Definition 2 および証明中の "non-adjacent cells" → "cells not connected by an edge of $T$"（BFS 木エッジで非隣接，という意味を明示） |
| M1: 等変性 Remark の不正確な例を削除 | Remark 3 から "for example, lexicographic order on cell centroids after the face-up rotation" を削除（SO(4) 回転は辞書順を保存しないため等変でない） |
| M2: Remark 1 の移動 | validity test の注釈を Definition 1 (BFS spanning tree) 直後から Definition 2 (Valid net) 直前に移動（論理的文脈の一致） |
| M2': h_i 記法の定義追加 | Mathematica コード後に "where $h_i$ denotes the convex hull mesh of cell $i$'s vertices in 3D" を追加（実装段落で未定義の記号 $h_i$ を使用していたため） |
| M3: Introduction 節参照の整理 | subsection 参照（sec:background, sec:algorithm, sec:contrast）を除去し `sec:theory`（§2 全体）と `sec:results`（§3）の上位節参照に統一 |
| L1: Introduction の重複定理を削除 | Theorem 1（Introduction）を削除し非形式的1文に置換；§3.1 の Theorem "Main, restated" を "Main"（ラベル `thm:main`）に改名；旧 `thm:main2` は消滅 |
| L2: コメント残骸の削除 | `% \paragraph{Relation to apple-peel unfolding.}` を削除 |
| L3: 括弧書きの修正 | "(the ordering always succeeds)" → "(completing for all 1,440 starting pairs)"（"valid" の二義性を排除） |
| L4: 証明中の "non-adjacent cells in BFS net" 修正 | → "cells not connected by a tree edge"（H2 と整合させる） |
| L5: "valid orderings" の二義性解消 | Introduction 中の "1,440 valid orderings (completing for all 1,440 starting pairs) but all resulting 3D nets self-intersect." → "That algorithm completes successfully for all 1,440 starting pairs for the 120-cell, yet every resulting 3D net self-intersects."（"valid" が ordering の完了と net の非自己交差の両方に使われていた問題を解消） |
| L6: Table 2 キャプション修正 | "(root, second cell) pairs" → "$(C_1, C_2)$ pairs"（記法を本文と統一） |
| L7: abstract の "diverges from algorithm" 修正 | "the result diverges sharply from the apple-peel unfolding algorithm" → "our result stands in sharp contrast to the apple-peel unfolding algorithm"（より正確な対比表現に） |
| L8: 未使用の face set $F$ を削除 | "vertex set $V$, face set $F$, cell set $C$" → "vertex set $V$, cell set $C$"（$F$ は Preliminaries 節以降で使用されていないため） |
| L9: "from the standard tables" 修正 | → "from the data files available in the repository"（Coxeter (1973) に 120-cell 頂点座標のミスプリントがあるため，使用したデータファイル `_4DData/f120.m` 等が修正済みであることを明示） |

#### 2026-08-19 セッションの修正（runcinated 24-cell の訂正 — **JoCG 投稿版**）

**方針（ユーザー判断 2026-08-19）**：対称性の理論は入れず，この論文はここまででまとめて **JoCG に投稿**する。対称性は次の論文へ。

| 修正 | 内容 |
|------|------|
| **§4 に段落追加（必須の訂正）** | *A net where no rule of Table 3 finds one.* — runcinated 24-cell は index 系3規則すべてで 0/240 だが，**等変な木（\|H\|=3）なら全240根から valid なネットが存在する**。Aut = 2304（拡大 [3,4,3]，図が回文），セル軌道2つ（48/192），各軌道の根で600本サンプルして 70本・83本 valid。群が各軌道に推移的なので**全根に伝播**。対照：無制約ランダム600本×2根で valid ゼロ |
| §4 の既存記述に前方参照 | 「yields one from no root」の直後に「これは $\mathcal{T}_r$ についての言明であって多胞体の性質ではない，節末を見よ」を挿入。**Open Problems の問いが誤った前提に立っていたのが訂正の動機** |
| Abstract に1文 | grand antiprism より**鋭い例**として runcinated 24-cell を追記 |
| §5 Open Problems を書き換え | 「不変量の候補が無い」→ **何が効かないかを名指し**：セル次数に加えて**辺の角度欠損**を排除。正多胞体6種を難しさ順に並べる（148.4/90/77.9/31.6/10.3/7.4°）のに一様では両方向に反例（x5x3x3x が MIXED と同値 0.736° で全根 valid，x3o4o3x は 10.5° で invalid）。深さ・サイズで重み付けても直らず正多胞体の順序が壊れることも記載 |
| §3 に Remark 追加 | *The 8-cell has only one net.* 層が 1+6+1 に固定 → 6本の木は6向きの **Dalí の十字架**で，48ネットすべて合同 |

11 → **12ページ**でビルド通過。未定義参照・引用警告なし。Overfull/Underfull 4件はすべて既存（Table 1 の幅と Code and Data 節）。

#### 2026-08-23 セッションの修正（なぜ validity を要求するのかの動機付け）

**発端**：ユーザーの問い「なぜ valid にこだわるのかの説明が論文にあったか」。**無かった**（§2.1 で "a net if no two cells overlap in their interiors" と定義するだけ。`paper_draft.tex` も Conclusion に "3D-printability" が一度出るのみ）。

| 修正 | 内容 |
|------|------|
| §2.1 に Remark 追加（`rem:whyvalid`） | *Why overlap must be excluded.* net の定義直後，Stella4D 段落の前。3段構成：(1) 順方向の写像は重なりの有無によらず定義でき逆回転で $P$ に戻る，(2) 壊れるのは**逆の割り当て**で，重なり点は2胞に属し折ると $P$ の相異なる2胞へ送られる → 1個の物質が $\mathbb{R}^4$ の2箇所を同時に占めることになり不可能，(3) 切り出す側から言えば重なりのある配置は $H$ の部分集合ではなく**はめ込みの像**なので1個の立体ブロックから切り出せない（3D で1枚の紙から切り出せることを要求するのと同じ） |
| Introduction に1文 | 「4次元では胞を $\RR^3$ に並べたものがネット」の直後に「重なりの排除は見た目の問題ではない：2回覆われた点は折り上げると2つの異なる胞へ行かねばならず，単射な展開だけが1個の立体から切り出せる」を1文。詳細は `Remark~\ref{rem:whyvalid}` へ前方参照 |

**12ページのまま**ビルド通過。未定義参照・引用警告なし。Overfull/Underfull は既存4件のみで増減なし。

> ### ⚠️ 「valid でないと元の多胞体に戻らない」は**誤り**
>
> 当初ユーザーが提案した理由付けはこれだったが，展開は木の各辺ごとの剛体回転の合成なので，**逆回転を辿れば重なりの有無によらず必ず $P$ に戻る**（3D でも，自己重なりのある展開図は数学的には折れば元の多面体になる）。この理由で書くと査読で確実に突かれる。
>
> 正しくは「**戻せるか**」ではなく「**展開写像 $\partial P \to H$ が単射か**」＝埋め込みか単なるはめ込みか。ユーザーの直観のうち救えるのは「**点集合としてのネットからは $P$ を復元できない**」という言い方（木のラベル情報を使ってよいなら常に戻る）。
>
> **紙の厚みについての留保**（本文には書いていない）：現実の紙は重ねて折れてしまうが，それは「厚みぶんだけ $w$ 方向にずらす」ことに相当し，$H$ 内の厚さゼロの立体という前提を捨てている。
>
> **覚え方（ユーザー，2026-08-23）**：「重なった紙の図形を，辺で折れば立体になるからといって**展開図と呼ぶか**？」— 折れるかどうかは初めから争点ではなく，1枚から切り出せる型紙かどうかが争点。

#### 2026-09-06 セッションの修正（内部矛盾3件 — **12 → 13ページ**）

発表準備でスライドと論文を突き合わせて見つかったもの。いずれも**論文内部の矛盾**。

| 修正 | 内容 |
|------|------|
| **(A) `rem:equivariance` の最終文** | "Whether validity holds for every possible tie-breaking rule remains open." は part 1 が 5/8/16/24胞体について**決着させている**ので誤り。**120/600胞体に限定**。あわせて論理を整理：**part 1 に等変性の仮定は要らない**（対称性が根 $r$ の BFS 木の**族全体**を根 $r'$ の族へ全単射で写すので，規則が何であれ1根で足りる）。等変性の但し書きが要るのは1根1本の $\mathcal{T}_r$ を名指しする **part 2 だけ** |
| **(B) §4 の random-tree の文** | "as does the grand antiprism, even though each of them is valid from every root under both index rules" は同じ論文の Table 3 と矛盾（gap は **min-index で 0/320**）。さらに 2,640セルの4種のうち**2種は MIXED なので「61種」の外**。両方を書き分けた |
| **(C) Seamons 2026 の反映** | §3.1 に `Remark 7`（`\label{rem:notallnet}`）を新設，§5 の第2 open problem を書き換え（**24胞体は否定的に決着**），Abstract に1文，`\bibitem{Seamons2026}`（Zenodo, doi:10.5281/zenodo.22004769）追加で **References 7 → 8件** |

Abstract の "By contrast, our result stands in sharp contrast to..." の重複も1語削った。
未定義参照なし，Overfull/Underfull は既存4件のみ。

> **2,640セルの多胞体はちょうど4種**（`_4DData_uniform/index.csv` で確認）：runcinated 120-cell `x5o3o3x`／
> runcitruncated 600-cell `x5o3x3x`／runcitruncated 120-cell `x5x3o3x`／omnitruncated 120-cell `x5x3x3x`。
> **真ん中の2つが MIXED**。「2640セル級」とまとめて書くと61種の内訳と食い違う。

> **Remark 7 の中身**：Seamons の木が BFS でないことの根拠を数字で示した。24胞体は任意の根から
> 離心数3・層 $1+8+14+1$ なのに，彼の木では **13/24 セルが距離より深い**。「preprint での announce に
> 留まる／著者が独立に厳密算術で再現した」の留保も入れた。八面体セルが中心対称ゆえ各蝶番が共有三角形の
> 平面での鏡映になり，展開が有理数に留まる。

#### 2026-08-11 セッションの修正（Stella4D の引用追加）

| 修正 | 内容 |
|------|------|
| §2.1 に段落追加 | ridge unfolding / net の定義直後（`\begin{definition}[BFS spanning tree]` の前）に Stella4D への言及を 1 段落。「既存ソフトはセル間の交差を無視するため，自己交差の有無は射程外」という対比 |
| 参考文献追加 | `\bibitem{Webb2026}`（R. Webb, *Stella4D: Polyhedron Navigator*, version 6.0, Software3D, 2026）を Turney と YoshinoChaidee の間に追加 |
| §4 Open Problems の最終段落を差し替え | 漠然とした "non-regular convex 4-polytopes (e.g. the regular-faced uniform 4-polytopes)" → **凸一様多胞体64種**への拡張を具体的に記述．(1) 無限系列2つを除き丁度64種で**完全性が証明済み**（星型込みの全一様枚挙は未解決なのと対比），例外的な非Wythoff的メンバーが grand antiprism，(2) 3D の正多面体→アルキメデス立体拡張（companion paper）の4次元版にあたる，(3) 障害は概念的でなく計算量的：最大の omnitruncated 120-cell はセル数 **2640** で全根の交差判定は高コスト，単一根またはサンプリングで妥協が要る |
| 参考文献追加 | `\bibitem{ConwayGuy1965}`（Conway & Guy, *Four-dimensional Archimedean polytopes*, Proc. Colloquium on Convexity, Copenhagen, 1965, pp. 38--39）を Coxeter と DevadossHarvey の間に追加．References 5件 → **7件** |

pdflatex 2 回通し済み。未定義参照・引用警告なし。**6ページ → 7ページ**になったが，7ページ目は参考文献の末尾3件のみで本文は6ページに収まっている。ドラフト段階ではページ数にこだわらない方針（ユーザー判断，2026-08-11）。6ページに戻したい場合は Open Problems の段落を3〜4行削れば足りる（omnitruncated 120-cell の計算量の記述が最も落としやすい）。

##### 一様多胞体の数（背景知識）

| | 正則 | 凸一様 | 一様（星型込み） |
|---|---:|---:|---:|
| 3D | 5 | 13（アルキメデス）+ 角柱・反角柱 | **75**（+無限系列）— 完全性証明済み（Sopov 1970 / Skilling 1975） |
| 4D | 6 | **64**（+無限系列2つ）— 証明済み（Conway & Guy 1965 の grand antiprism が例外的非Wythoff） | **2191（既知）— 枚挙は未解決** |

2191 は「証明された総数」ではなく既知数。枚挙はほぼアマチュア主導（Jonathan Bowers，George Olshevsky）で，2021年1月に新 snub regiment が272個見つかり 2127，同年4月に333個へ育って 2188，最後の2個が2023年4月。Stella 5.4 (2014) の 1849 → 6.0 の 2191（+342）はこの12年の発見を反映。changelog の *fissary* / *regiment*・*coincidic* / *scaliform* は Bowers の造語で，査読論文の標準語彙ではないため使用時は注意。

#### Stella4D（4D 多胞体の 3D ネット表示ソフト）

Robert Webb 作，<https://www.software3d.com/Stella.php>。4D 多胞体の 3D ネットを表示できる唯一の主要ソフト。
**本論文の貢献は先取りされていない**（マニュアル §15.8 "4D Nets" の記述による）：

- "attempts to generate nets with as much symmetry and **aesthetic appeal** as possible" — 展開木の族（BFS 等）を指定する仕組みではない
- "in 4D ... **intersections between cells are ignored**" — **自己交差・重なりの判定機能がない**
- `Ctrl+右クリック` で隣接セルを貼り替え，`Nets > Maximum Cells per Net = 1` で 1 セルずつ育てられる（任意 spanning tree を手動でなら作れるが自動列挙はしない）
- §15.5：4D では "The symmetry group of a 4D polytope is **not established**"（ヒューリスティック分類のみ）

**Stella 6.0**（2026-08-11 リリース，12年ぶりのメジャー更新）の changelog にも 4D ネットの重なり判定・列挙に関する項目は皆無（ネット関連の追加はすべて 2D ネット向け）。有用な変更点：

- **64-bit 化。Mac/Linux の Wine で良好に動作**（Parallels 不要）。ただし M2 / Sonoma 以降では Wine のバグで印刷不可（表示・エクスポートは可）
- **4D OFF ファイルのインポートがスケールを保持するよう修正** → `_4DData/f120.m` 等を OFF 出力して読み込ませ，Stella4D の既定ネットと自作 BFS ネットを比較するクロスチェックが現実的に（**未実施**）
- 一様多胞体ライブラリ 1849 → 2191 種；`4D > Create Segmentotope`；投影方向 "Cell Last" 追加
- 価格：新規 US$67（Pro $120），既存ユーザーの 6.0 アップグレード US$20

**方針（2026-08-11 決定）**：**ソフトは購入しない**。手元にないのは一様多胞体ライブラリ（論文は正凸 6 種のみ）と既定ネットの見た目だけで，しかも Stella4D は**どの spanning tree を使ったかを表示しない**ため買っても「既定ネットは BFS 族か」に即答できない。論文の主張は Stella の既定ネットに依存しない。
**査読者が Stella4D との比較を求めてきた場合に限り**，作者 Robert Webb にサイトの Contact 経由で「既定の 4D ネットの生成規則は何か，BFS 的なものか」を問い合わせる（無料・権威ある回答・*personal communication* として引用可）。購入検討はそれでも足りない場合のみ。

> **注意**：software3d.com は bot の User-Agent に HTTP 403 を返す。WebFetch は失敗するので `curl -A "Mozilla/5.0 ..."` で取得すること。

#### Coxeter (1973) ミスプリントについて

Coxeter "Regular Polytopes" 3rd ed. (Dover, 1973) に 120-cell の**頂点座標**のミスプリントが存在する。組合せデータ（600頂点・1200辺・720面・120セル）は正しい。本論文で使用したデータファイル（`_4DData/f120.m` 等）はミスプリントを修正済みのため，計算結果への影響はない。これが L9 で "from the data files available in the repository" とした理由。

#### arXiv 投稿状況（2026-07-29 更新）

- zip ファイル（`face_rotation_net_paper.tex` + `260612AllFaceRotation.pdf`，計 249KB）を作成し arXiv に投稿
- **submission ID：submit/7751072**（primary category = **cs.CG**）
- 2026-06-25 投稿 → **on hold のまま34日以上経過**（ステータスは純粋に "on hold" のみ。endorsement 系ではなく moderation で停止）
- companion 論文 **2605.30373 も cs.CG で公開済み**（同じ著者・同じカテゴリが通っている → 「cs.CG だから遅い」ではなく，この個別 submission が moderator のところで放置されている）
- help@arxiv.org に問い合わせ済み。**2026-07-15 に返信あり**だがテンプレート（「moderator に催促した，対応不要」）のみ。以後も動かず
- **再問い合わせメール：2026-08-07 送信済み**（JCDCG³ 採択を新しい判断材料として追加した文案。submit/7751072・2604.16204・2605.30373 を明記）。**受け取り確認の返信はあったが，2026-08-11 時点で on hold のまま変化なし**
- 方針：この1通だけ送ってあとは完全放置。**arXiv プレプリントは JCDCG³ 発表・journal 投稿の前提ではない**ため，hold が続いても本筋の損失はない
- **JCDCG³ 2026 は採択済み**（2026-07-03 通知）→ proceedings 掲載・発表は確保。journal 投稿（第1候補 JoCG）は arXiv の状態と無関係に進められる

#### "We" vs "I"（単著慣習）

単著論文で "We" を使用（11箇所）。数学・計算幾何分野では単著でも authorial "we" が標準慣習のため，そのまま維持することにした。

---

## 凸一様4-多胞体への拡張（Wythoff 構成，2026-08-12）

face-rotation BFS ネットの検証を，正多胞体6種から**凸一様多胞体**へ拡張。基本データを Wythoff 構成で自前生成し，既存の BFS パイプラインに流した。

### スクリプト

| ファイル | 役割 |
|----------|------|
| `src/uniform/uniform4D_wythoff.m` | **生成器**：Coxeter 図から一様4-多胞体の頂点・辺・面・セルを生成（A₄/B₄/F₄，33種） |
| `src/uniform/run_wythoff_bfs.m` | 生成データを BFS チェッカに流す（Phase A：正多胞体を全根／Phase B：全種を root 1） |
| `src/uniform/run_wythoff_bfs_allroots.m` | **全30種 × 全根**（1893根）の確定計算。結果を `data/wythoff_bfs_allroots.mx` に保存 |
| `check_runcinated24.m` | 唯一の例外 x3o4o3x の全240根詳細（重なりペア数の分布） |

`src/bfsnet/face_rotation_net_all4D_v2.m` は本体の main ループを走らせずに定義だけ読み込む（テキストを `Print["Face-rotation BFS net check v2` の直前で切って `ToExpression`）。

### 生成アルゴリズム（`src/uniform/uniform4D_wythoff.m`）

1. Coxeter マーク → Gram 行列 → Cholesky で単位鏡法線 n_i
2. 群 W = 4つの鏡映の閉包（ハッシュキーは `Round[x/grid]` で**厳密整数**にする。`Round[x, dx]` は精度付き実数を返しキーに使えない）
3. 種点 p は `n_i . p = b_i`（b_i = 1 リング付き／0 なし）。**鏡 i が生む辺の長さは 2b_i** なので，b_i を揃えれば辺長が自動的に揃う＝一様性が担保される
4. 頂点 = p の軌道
5. **セル**：ノード i を除いた極大放物型部分群 W_J が facet を固定する。その超平面を W で軌道に乗せ，各超平面上の頂点を集める（凸包計算は不要，厳密かつ高速）
6. **2-面（ridge）**：凸4-多胞体の ridge はちょうど2つの facet の交わり → セル対の共有頂点集合でアフィン次元2のもの
7. 面は**巡回順**で格納（`unfoldTFP` が共有面の最初の2頂点から回転平面を作るため，正方形を添字順に並べると対蹠点になり退化する）

### 自己検証（文献値に依存しない）

- **V** = |W| / |W_unringed|（種点の固定部分群はリングなしノードの放物型部分群）
- **C** = Σᵢ |W| / |W_Jᵢ|（facet を持つノードについて和）
- **オイラー標数** χ = V − E + F − C = 0
- 全 facet 超平面が実際に多胞体を支持しているか

**33/33 が全4項目を通過**（群位数も 120 / 384 / 1152 と既知値どおり）。

### 生成結果

| 族 | Wythoff 的メンバー | |W| |
|----|:---:|---:|
| A₄ | 9 | 120 |
| B₄ | 15 | 384 |
| F₄ | 9 | 1152 |

族をまたぐ重複3組を自動検出（組合せシグネチャ一致）：`o4o3x3o == o3o4o3x`（24胞体），`o4x3o3x == o3o4x3o`，`o4x3x3x == o3o4x3x`。33 − 3 = **正味30種**。

**未対応**：snub 24-cell（交替構成）と grand antiprism（600胞体から直交2十角形リング計20頂点を除いた凸包）は Wythoff 的でないため別構成が必要。H₄（位数14400）は marks 差し替えで動くはずだが，ridge 検出が O(C²) のため omnitruncated 120-cell（2640セル）には書き換えが必要（rank-2 放物型から ridge を作り，重心×法線の行列積で所属セルを引く方針をヘッダコメントに記載）。

### BFS ネット計算結果（全30種 × 全根，1893根）

```
valid roots        : 1653 / 1893
ALL VALID polytopes: 29
MIXED polytopes    : 0        ← 中間が存在しない
ALL INVALID        : 1
  x3o4o3x   0/240   cells={{6,240}}   deg={{5,192},{8,48}}
```

- **29種は全根で valid，重なりペアは min/max とも 0**
- **x3o4o3x（runcinated 24-cell）だけが全240根で invalid**。最良の根でも重なり3組，最悪23組（分布 `{{3,7},{4,38},{5,36},{6,32},{7,7},{13,24},{21,51},{22,42},{23,3}}`）
- この30種の範囲では MIXED がゼロ。**ただしこれは H₄ で破れる**（下の「H₄ 系」参照）。当時「完全な二分法」と記録したが**反証済み**

### 重要な観察

例外の特徴：同じ240セルでも x3o4x3x・x3x4x3x は全根 valid。x3o4o3x は**最小セル次数5，かつ次数5のセルが8割**。600胞体（次数4）が apple-peel で impossible だった件と方向性は一致するが，次数5を含む o3x3o3x 等も valid なので**次数だけでは説明できない**。snub 24-cell と grand antiprism（最小次数4）も全根 valid なので，次数の低さは自己交差の原因ではない。

### Phase A：データ経路の検証

生成データで正多胞体を回し，既存 `_4DData/` と同一の結果を再現：5-cell 5/5，8-cell 8/8，16-cell 16/16，24-cell 24/24 いずれも ALL VALID。生成器 → BFS → 体積判定の経路全体が正しい。

### `_4DData/f*.m` の形式（当初の誤読を訂正）

| 位置 | 内容 |
|:----:|------|
| [[1]] | 頂点座標 |
| [[2]] | 辺（頂点添字ペア） |
| [[3]] | **各頂点の隣接頂点リスト**（セルの頂点リストではない） |
| [[4]] | 2-面（頂点添字リスト） |
| [[5]] | 面隣接（辺を共有する面） |
| [[6]] | セル（面添字リスト） |

`src/bfsnet/face_rotation_net_all4D_v2.m` が使うのは [[1]], [[4]], [[6]] のみ。

- [[3]] を「セルの頂点リスト」と誤読しやすい：f5.m では頂点数5・次数4がセル数5・4頂点と**偶然一致**する
- **`f8.m` の [[2]] は32本の辺を両方向で64件**列挙（f5/f16/f24 は1回ずつ）。[[2]] は誰も読まないので実害なし。生成側は無向1回の規約

### H₄ 系（2026-08-12 追加）

ridge 検出を **O(C²) → O(|W|·多角形サイズ)** に書き換えて H₄（位数14400）に到達。

**書き換えの中身**（`src/uniform/uniform4D_wythoff.m`）
- **ridge**：Wythoff 的多胞体の d-面は階数 d の放物型部分群の軌道。rank-2 放物型 W_{i,j} で ridge の原型を作り，その**頂点集合**を W で軌道に乗せる。セル対の総当たりが不要に
- **ridge が属する2セル**：頂点→セルの逆引きを作り，ridge の全頂点が属するセルの共通部分（ちょうど2個になるはず）
- **セル所属判定**：超平面ごとに1回の行列ベクトル積にベクトル化
- 検証項目が2つ増加：**各 ridge がちょうど2セルに属する**，**全頂点の次数が一様（頂点推移性）**
- A₄/B₄/F₄ の33種は書き換え前と**完全に同一の V/E/F/C**（回帰テスト通過）かつ高速化

**生成結果：48/48 が全自己検証を通過**（A₄ 9 + B₄ 15 + F₄ 9 + H₄ 15），重複3組を除いて**正味45種**。総生成時間509秒，`data/wythoff_gen_all.mx` に保存。

| symbol | V | E | F | C | |
|---|---:|---:|---:|---:|---|
| x5o3o3o | 600 | 1200 | 720 | 120 | 120胞体 |
| o5o3o3x | 120 | 720 | 1200 | 600 | 600胞体 |
| x5x3x3x | 14400 | 28800 | 17040 | 2640 | omnitruncated 120-cell |

`f120.m` / `f600.m` との照合も MATCH（辺数含む）。

**BFS 全根（2026-08-13，SAT 版で完走）：21,600根**

```
polytopes 15   roots 21600   valid roots 18937
ALL VALID 13   MIXED 2   ALL INVALID 0

x5o3x3x   2640 cells     47/2640   MIXED   overlaps 0/11
x5x3o3x   2640 cells   2570/2640   MIXED   overlaps 0/1
```

> ### ⚠️ 二分法は反証された（重要な訂正）
>
> **A₄/B₄/F₄ の30種と snub/grand antiprism で観察された「全か無か」は H₄ で破れる。** 上記2種が MIXED。「MIXED ゼロ＝完全な二分法」という以前の記述は**誤りなので参照しないこと**。
>
> **「root 1 は強い指標」も誤り。** x5o3x3x は root 1 で VALID だが全根では **47/2640** しか valid でない（root 1 がたまたま47個のうちの1つだった）。root 1 のみの結果に判断材料としての価値はほとんどない。

**MIXED の裏取り（`src/uniform/verify_mixed_h4.m`，2026-08-13）**：数値誤差の疑いを排除済み。

- 失敗根5個・成功根3個 × 2多胞体 = 16根で **SAT と RegionMeasure が判定・重なり数とも完全一致**
- eps を 10⁻⁴ 〜 10⁻¹⁰ まで振っても重なり数は不変（閾値上の縁の判定ではない）
- x5o3x3x の重なり数分布は滑らか：`{{0,47},{1,282},{2,361},{3,546},{4,608},{5,439},{6,225},{7,76},{8,34},{9,12},{10,9},{11,1}}`
- x5x3o3x は `{{0,2570},{1,70}}` で，失敗根はすべて重なり1個ちょうど

### snub 24-cell と grand antiprism（2026-08-12，`src/uniform/uniform4D_special.m`）

Wythoff 構成（リング付き Coxeter 図）で到達できない2種。**どちらも「600胞体から頂点を除いて凸包を取る」同一の構成**で作れる。

- 除いた各頂点 v について，v の12個の隣接頂点のうち**生き残った分**が新しいセルになる
- **snub 24-cell**：除くのは内接24胞体（24頂点，600胞体内で互いに非隣接）→ 各キャップは12頂点全部残って**二十面体**
- **grand antiprism**：除くのは直交する2つの十角形リング（20頂点）→ 各キャップはリング上の隣2個が消えて10頂点＝**五角反柱**

600胞体は Wythoff 生成器からではなく**標準 icosian 座標**で直接構成する（内接24胞体が単に最初の24頂点＝±e_i と (±½,±½,±½,±½) になり，探索が不要になるため）。セルは辺グラフの**4-クリーク**（各頂点のリンクが二十面体で面20個 → 120×20/4 = 600）。十角形リングは72本見つかり既知の値と一致。

```
polytope          symbol      V     E     F    C   status
snub 24-cell      s3s4o3o    96   432   480  144   OK  cells {{4,120},{12,24}}  deg {{9,96}}
grand antiprism   gap       100   500   720  320   OK  cells {{4,300},{10,20}}  deg {{10,100}}
```

**BFS 全根：両方とも ALL VALID，重なり0**（144/144 が31.6秒，320/320 が344.9秒）。結果は `data/special_bfs_allroots.mx`。

### 角柱17種（2026-08-13，`src/uniform/uniform4D_prisms.m`）

Platonic 5 + Archimedean 13 − 立方体角柱（＝8胞体，正則なので既出）= **17種**。3D 多面体 P の角柱は**組合せ的に直接構成でき，凸包計算は不要**：

- 頂点 2n（V(P) の2コピー），セル F+2（P の2コピー + 各面ごとの角柱），2-面 2F+E，辺 2E+n
- オイラー：2n−(2E+n)+(2F+E)−(F+2) = (n−E+F)−2 = 0 ✓
- 一様性のため側面の辺長 = P の辺長（P を辺長1に正規化して高さ1）
- 3D データは `PolyhedronData`。incidence は `sp4Incidence`（C ≤ 94 なので O(C²) で十分）

**17/17 が検証を通過。BFS 全530根すべて ALL VALID，重なり0。**

角柱のセル隣接グラフは2つの「蓋」セルが全側面セルに隣接するハブ構造で，蓋の次数は最大92（snub dodecahedron 角柱）。それでも全根 valid。

### 凸一様4-多胞体の全体像（2026-08-13，**全64種を全根で確定**）

| 対象 | 種数 | 根数 | 結果 |
|---|:---:|---:|---|
| A₄/B₄/F₄（正味30種） | 30 | 1,893 | 29種 ALL VALID／x3o4o3x のみ 0/240 |
| H₄ | 15 | 21,600 | 13種 ALL VALID／x5o3x3x 47/2640・x5x3o3x 2570/2640 が MIXED |
| snub 24-cell + grand antiprism | 2 | 464 | 2種とも ALL VALID |
| 角柱 | 17 | 530 | 17種とも ALL VALID |

**凸一様64種・全24,487根が確定**（valid 21,584根）。

| 分類 | 種数 | 内訳 |
|---|:---:|---|
| ALL VALID | 61 | |
| **MIXED** | **2** | x5o3x3x = $t_{0,2,3}\{5,3,3\}$ runcitruncated 600-cell (47/2640)<br>x5x3o3x = $t_{0,1,3}\{5,3,3\}$ runcitruncated 120-cell (2570/2640) |
| ALL INVALID | 1 | x3o4o3x = $t_{0,3}\{3,4,3\}$ runcinated 24-cell (0/240) |

> **命名の注意**：x5o3x3x と x5x3o3x は取り違えやすい。セル構成で判別すること。
> x5o3x3x = 120 truncated icosahedra + 720 pentagonal prisms + 1200 hexagonal prisms + 600 cuboctahedra（**600-cell** 側）。
> x5x3o3x = 120 truncated dodecahedra + 720 decagonal prisms + 1200 triangular prisms + 600 cuboctahedra（**120-cell** 側）。

**セル次数は両方向に効かない**：低次数でもダメではない（snub 24-cell・grand antiprism は最小次数4だが全根 valid），高次数でも保証にならない（角柱は次数92のセルを持ち全根 valid だが，MIXED の x5x3o3x も次数32のセルを持つ）。x3o4o3x は最小次数5で唯一の全根 invalid。**3つの結果を分ける不変量は不明。**

### 組合せデータ公開（`_4DData_uniform/`，2026-08-15）

`src/uniform/export_uniform_data.m` で64種のデータを出力。**64ファイル 13.3 MB**。`index.csv` 付き。

**5要素形式**（`_4DData/f*.m` の6要素から2つ削り1つ足した）：

```
[[1]] verts        頂点座標（機械精度）
[[2]] edges        頂点添字ペア
[[3]] faces        2-面（頂点添字，巡回順）
[[4]] cellsByFace  セル（面添字）
[[5]] cellAdj      各セルの隣接セル  ← 新規
```

**削った2つ**：`vertAdj`（第3要素）と `faceAdj`（第5要素）。理由は容量ではない。
- **`src/bfsnet/face_rotation_net_all4D_v2.m` は `raw[[1]], raw[[4]], raw[[6]]` の3つしか読まない**（138行目）。この2つはリポジトリ内のどこからも参照されていない
- 45種で 2.2 MB と 6.4 MB。README に**1行での復元コード**を掲載済み

**足した1つ**：`cellAdj`。**ネットに必要な唯一の隣接情報**であり，`_4DData/` はこれを持っていない（パイプラインは読み込みのたびに `cellsF` から再計算している）。45種で 1.3 MB＝面隣接の5分の1。

**座標は機械精度**（生成器の30桁ではない）。パイプラインが読み込み直後に `N[...]` で落とすため，30桁は誰も使っておらず，配布すると「この精度に意味がある」と誤解させる。

書き出し後に**読み直して検証**：オイラー標数，セル隣接の対称性，セル隣接が「面添字を共有」の定義と厳密一致。

> **サイズについての整理（2026-08-15）**：12MB 程度では git / GitHub の仕様上の問題は**ない**。
> 実測で `uniform_nets/`（作業ツリー 56.7MB）は git パック内で **17.0MB（30%）**。GitHub の制限は
> 1ファイル100MB・リポジトリ推奨1GB で桁が違う。気にすべきは容量ではなく，(1) git は消せないので
> **生成物を作り直すたびに履歴へ永久追加される**こと，(2) clone コストが読者に転嫁されること。
> 上の「削る/足す」判断は**容量ではなく正しさと有用性**が理由。

### STL / OFF 公開（`uniform_nets/`，2026-08-13）

`src/uniform/export_uniform_nets.m` で **67ケース**（64種の root 1 ＋ 例外3種の追加根）を出力。合計55MB。

```
uniform_nets/
  stl/   67 files, 32 MB   バイナリSTL（セルごとに閉じたシェル，fan 三角形分割）
  off/   67 files, 23 MB   OFF（多角形面を保持，頂点は位置で重複除去）
  manifest.csv, README.md
```

- **OFF でセル構造を保存**：OFF にセルの概念はないので，面リストの後に `# cell k: <0-based face indices>` のコメント行で記録。標準ビューアは無視するがパースすれば復元可能
- **manifest の name 列は Schläfli t 記法**（リング配置から機械生成）。45種に英語名を付けると転記ミスが入るため
- **自己交差ファイルは4件**（x3o4o3x root1/root47，x5o3x3x root14，x5x3o3x root748）。**非多様体でスライサに通せない**旨を README で警告
- `.gitignore` は `*.stl` を全体除外しているので `!uniform_nets/stl/*.stl` 等の例外規則を追加済み
- STL の書式検証：ファイルサイズ = 84 + 50×三角形数 と厳密一致

**セル次数は自己交差の原因ではない**：snub 24-cell も grand antiprism も最小セル次数4（600胞体が apple-peel で impossible だったのと同じ「貧しさ」）だが face-rotation BFS では全根 valid。x3o4o3x の最小次数5がそれでも唯一の反例である理由は依然として不明。

### 論文改稿（2026-08-13，`face_rotation_net_paper.tex`）

**タイトル変更**：`Face-Rotation BFS Nets of the Six Regular Convex 4-Polytopes`
→ **`Face-Rotation BFS Nets of the Regular and Uniform Convex 4-Polytopes`**

**構成**（一様の結果を独立節に昇格）：

```
1. Introduction
2. Preliminaries   2.1 Background / 2.2 Face-Rotation BFS Algorithm
3. Results         3.1 Main result / 3.2 Contrast with Apple-Peel Unfolding
4. The Uniform Case  ← 新設（\label{sec:uniform}）
5. Open Problems
```

§4 は3つの paragraph：*Construction*（Wythoff 45種・600胞体の diminishing 2種・角柱17種の作り方と検証4項目）／*Results*（表 `tab:uniform` = 4群×5列，三分の内訳，根依存2種の議論，セル次数が両方向に効かない段落）／*Caveats*（タイブレーク規則依存の理由，SAT の妥当性検証）。

**Figure 2（`fig:uniformfail`）新設**：`src/figures/make_uniform_figure.m` で生成。失敗の2つの型を1つずつ。
- (a) runcinated 24-cell の**最良の根**（root 47，交差3組）を**全体図**で — 「どの根でもダメ」は最悪例より最良例が効く
- (b) runcitruncated 120-cell の失敗根（root 748）の**交差ペア近傍8セルの拡大図**— 2640セル全体は団子になるため
- 交差ペアの2セルは**赤と青に塗り分け**（同色だと相貫が1個の多面体に見える）

**2枚の縮尺が違うのは意図的**（2026-08-17 にキャプションで明示）。交差の局所的な見た目はどちらの多胞体でも同じなので，両方を拡大図にすると**同じ絵が2枚並ぶ**ことになり，(a) 0/240 と (b) 2570/2640 という違いが伝わらない。(a) を全体図にすると「240セルの中に交差が3組だけ散在し，それが240通りの根の最良の結果」という意外性が出る。**縮尺差そのものが主張の粒度（(a)=多胞体について，(b)=1つの根について）に対応している。** キャプションに "drawn whole" と "drawn as a close-up" を対比させて明示済み。

#### 図のスタイル（2026-08-15 確定，Figure 1 に合わせて黒背景化）

**Figure 1（`260612AllFaceRotation.pdf`）の実測**：背景は**純黒ではなく GrayLevel[0.13]**（濃いチャコール），セルは**不透明**，**エッジ線は描かれていない**（面の境界は照明の陰影差のみ）。

当初 Figure 2 を白背景で作ったため**同一論文内で Figure 1 と不整合**だった。黒版に統一。

| 項目 | 設定 | 理由 |
|------|------|------|
| 背景 | `GrayLevel[0.13]` | Figure 1 の実測値 |
| 文脈セル | エッジなし，`GrayLevel[0.92, op]` | Figure 1 は陰影のみ。暗いエッジも試したが線が累積してワイヤーフレーム化し質感が離れる |
| 交差ペア | **エッジなし**，赤 `RGBColor[1,0.28,0.24,0.95]` / 青 `RGBColor[0.35,0.62,1,0.95]` | 陰影のみで面が分かれ Figure 1 と揃う |
| 不透明度 | **パネルごとに変える**：(a) 0.13，(b) 0.25 | 重なるセル数に反比例させる必要がある。(b) の36セルで最適な 0.25 は (a) の240セルでは下半分が白い塊に飽和する |

`src/figures/make_uniform_figure.m` は white / black 両方を出力（`fig_uniform_a.png` / `fig_uniform_a_black.png` 等）。論文は black 版を使用。

**パネル(b) の表示範囲（2026-08-17 修正）**：当初 `PlotRange` を明示していたが，選択半径（重心 5·cellR 以内）と表示範囲（3.2·cellR）が食い違っており，**セルが箱の面でスライスされて切断面が見えていた**。

試して却下した案：

| 案 | 結果 |
|----|------|
| 箱に完全に収まるセルだけ描く | ネットは連結しているので**穴が空いて破片が浮く**（36→7個）。ネットに見えない |
| 選択セル全部が収まるよう箱を広げる | 遠くまで伸びる1個が範囲を支配し，ペアが小さくなって黒地が余る |
| **`PlotRange` を指定せず自動**（採用） | 描くセルに範囲が密着。切断ゼロ・余白ゼロ |

採用形：**重心 3·cellR 以内のセル（8個）を選び，`PlotRange` は指定しない**。重心で選ぶので連結性は保たれる。パネル(a) は全セルを描き `PlotRange` も指定していないため元から問題なし。

> **未解決の限界（2点）**
> 1. パネル(a) の3組目の交差ペアは深部にあり，不透明度をどう調整しても明瞭にならない。図の役割は「最良の根でも3組が交差する」ことなのでキャプションの数字で足りるが，3組すべてを見せたい場合は視点変更か別カットが必要
> 2. パネル(b) はパネル(a) より約17%縦長で，横並びにすると高さが揃わない。**これは余白ではなくトリミングできない** — 一見空白に見える下部にも淡い半透明セルが実在する（コンテンツマスクで実測，両者とも既に密着）

Abstract に1文追加（64種24,487根，61/1/2）。Introduction の構成説明に §4 を追加。Code and Data Availability に角柱・SAT・`uniform_nets/` を追記。**9ページ**でビルド通過，警告なし。

### タイブレーク規則の徹底調査（2026-08-17〜19）

**発端**：タイブレークが index-base で幾何学的根拠がない。査読で必ず問われるとユーザーが判断し，拡張可能性の議論から始まった。結果として**論文の主張が大きく変わった**。

#### なぜ「幾何学的な規則に替える」では解決しないか

セル v の親を選ぶとき，**v の親候補が Stab(r, v)（根 r と v の両方を固定する対称性）の1つの軌道になる**なら，その群で不変な量はどれも候補間で同じ値になる。max w も回転角も Det も全てタイになる（apple-peel 側で xy クロス積という二次基準が必要だったのと同じ現象）。解は「良い規則を探す」ではなく「主張から規則を消す」こと。

> ### ⚠️ 「等変な規則は原理的に決められない」は正多胞体全体では成り立たない（2026-08-19 訂正）
>
> セル隣接グラフの Aut を全列挙して Stab(r, v) 軌道を実測した結果（`GraphAutomorphismGroup`，根 r=1）：
>
> | | 親候補が2個以上のセル | うち候補が単一軌道 |
> |---|---:|---:|
> | 5胞体 | 0（そもそも木が1本） | — |
> | 8胞体 | 1 | 1 |
> | 16胞体 | 11 | 11 |
> | 24胞体 | 7 | 7 |
> | **120胞体** | **71** | **21** |
> | **600胞体** | **183** | **99** |
>
> **120胞体では 71 個中 50 個，600胞体では 183 個中 84 個で候補が2つ以上の軌道に分裂する。** つまりこの2種に限れば**等変な規則が候補を区別できる余地がある**（そのような規則を探したわけではない）。「原理的に決められない」と書いていたのは誤りで，**5/8/16/24胞体でのみ検証された主張**である。
>
> 全6種で |Aut(セル隣接グラフ)| = 対称群の位数（120 / 384 / 384 / 1152 / 14400 / 14400）と一致するので，組合せ的自動同型＝幾何的対称性として読んでよい。
>
> 検証スクリプト：**`src/symmetry/check_tiebreak_orbits.wls`**（`wolframscript -file` で実行，全6種で約20分）。
>
> **落とし穴**：`Graph[edges]` の `VertexList` は出現順であり `Range[n]` ではない。`GraphAutomorphismGroup` が返す置換はこの添字に作用するので，頂点名でそのまま `PermutationReplace` すると**黙って別の頂点を動かす**。`Graph[Range[n], edges]` と明示すること。最初これで 16胞体・24胞体まで「軌道が分裂する」という誤った結果を得た。
>
> **`face_rotation_net_paper.tex` はこの過剰な主張を含んでいない**（Remark 4.1 は「実装の index タイブレークは等変でない」「全規則で成り立つかは未解決」と書くのみ）。修正が必要なのは本ファイルの記述だけ。

#### BFS 木の正しい定式化

**BFS の層構造（根からの距離）はタイブレークに依存しない不変量。** よって BFS 木とは「各セル $v$ について $N(v) \cap L_{d(v)-1}$ から親を選ぶこと」そのもので，木の総数は積

$$\prod_{v \neq r} |N(v) \cap L_{d(v)-1}|$$

| | セル | 根あたりの BFS 木 |
|---|---:|---:|
| 5胞体 | 5 | **1**（タイブレークの余地なし） |
| 8胞体 | 8 | 6 |
| 16胞体 | 16 | 20,736 |
| 24胞体 | 24 | 32,768 |
| 120胞体 | 120 | 1.76×10⁴⁴ |
| 600胞体 | 600 | 1.78×10⁶⁵ |

#### 結果1：正多胞体は規則に完全に鈍感（`src/tiebreak/tiebreak_regular.m`）

**5/8/16/24胞体は全 BFS 木（1,118,261本）が valid**。120/600胞体は根あたり20本の一様サンプル（14,400本）も全て valid。**総計 1,132,661 本で invalid ゼロ。**

**24胞体は D&H の未解決ケースで，BFS 族全体について決着した。** 定理を2部構成に書き換え（part 1 = 全 BFS 木，part 2 = $\mathcal{T}_r$）。

#### 結果2：探索が到達できる木は族の一部（`src/tiebreak/order_induced_trees.m`）

**Lemma**：BFS 木が探索で生成可能 ⟺ 親選択が全順序で誘導される ⟺ 制約有向グラフ $p(v) \to q$ が非巡回。

```
 8-cell        6 /      6   到達可能  100%
16-cell    4,896 / 20,736   到達可能   23.6%
24-cell    2,592 / 32,768   到達可能    7.9%
```

**24胞体で検証した 786,432 本のうち探索到達可能は 62,208 本（7.9%）だけ。** 残る 724,224 本はどんな BFS 実装も返さない木。定理 part 1 は素朴な読みの**12.6倍の範囲**を主張している。

> **障害は「候補を2つ共有する兄弟ペア」ではない**（それは十分条件にすぎない）。16胞体はそのペアが0組なのに4分の3が到達不可能 — $v_1$ が $a{\prec}b$，$v_2$ が $b{\prec}c$，$v_3$ が $c{\prec}a$ で**循環**するため。

#### 結果3：一様多胞体は規則に極度に敏感（`src/tiebreak/tiebreak_all61.m`, `src/tiebreak/tiebreak_maxindex.m`, `src/tiebreak/tiebreak_rules.m`）

```
                          queue        max-index     min-index
Runcinated 24-cell        0/240        0/240         0/240
Runcitruncated 600-cell   47/2640      2635/2640     2627/2640
Runcitruncated 120-cell   2570/2640    2640/2640     2151/2640
Grand antiprism           320/320      320/320       0/320     ← 決定的
```

**grand antiprism は queue と max-index で全320根 valid，min-index で1根も valid でない。** 「全根で valid」は多胞体の性質ではなく**(多胞体, 規則) の組の性質**。他に3種が根を失う（$t_{0,2,3}\{3,3,3\}$ が max-index で1根，$t_{0,3}\{4,3,3\}$ と $t_{0,3}\{5,3,3\}$ が min-index で1根・26根）。

**一様ランダム木では61種中13種が落ちる**。2640セルの多胞体は4種とも 0/10 だが，**うち2種（runcitruncated 600-cell と 120-cell）は MIXED なので61種には入っていない** — まとめて「61種のうち」と書かないこと（2026-09-06 に論文とスライドで修正）。grand antiprism も 0/10 で，こちらは queue と max-index では全320根 valid（min-index は 0/320 なので「両 index 規則で全根 valid」も誤り）。対して正多胞体は 1,132,661 本で失敗ゼロ。

> **結論：正則性とは「タイブレークに鈍感であること」だった。** 当初の三分類（61/2/1）より内容のある主張。

#### 私が犯した誤りの記録（再発防止）

| 誤り | 実際 |
|------|------|
| 「index 順がたまたま当たりを引いていた」 | 逆。queue 規則が**特に悪い**部類（x5o3x3x で 1.8%，min-index は 99.5%） |
| 「根依存性は根の性質」 | 誤り。**規則の性質**。rule-valid な根もランダム木ではほぼ全滅 |
| `MinimalBy[list, {f1,f2}]` で複数基準 | **この構文は存在しない**。未評価の式を比較し max-w/min-w が黙って index 順に退化していた |
| 「min-index は queue と同じ規則」 | 別物。queue の親は**最初にキューに入った候補**で，min-index とは 55%/25% 食い違う |
| Definition に「層条件 = 探索が生成できる木」と記載 | **偽**。層条件のほうが真に広い（24胞体で 92% が到達不可能） |
| 「gap 以外に2種が根を失う」 | **3種**。うち1種の名前も誤記（cantitruncated → runcitruncated 5-cell） |

**教訓**：30種・2,357根という規模でも，そこから外挿した性質（二分法）と近道（root 1 で代表させる）は両方とも外れた。

#### 論文への反映（すべて実装済み）

| 箇所 | 内容 |
|------|------|
| Definition 1 | 規則を定義から排除。BFS 木 = 最短路木。$G(P)$ が根に依存しないことも明記 |
| Lemma 1 | 探索到達可能性の特徴づけ＋証明 |
| Remark 1 | Table 2（到達可能な木の計数），循環の例 |
| Remark 2 | 実際に計算した木 $\mathcal{T}_r$ の特定と，それが素朴に見えない理由 |
| Theorem | 2部構成（part 1 = 全 BFS 木，part 2 = $\mathcal{T}_r$） |
| §4 | Table 3（3規則の比較），三分類が (多胞体, 規則) の分類であることを明示 |
| Table 1 キャプション | 同上の注記 |
| Abstract | 61/1/2 を撤回し**規則依存性**を主張の中心に |

> **記号注意**：木は $\mathcal{T}_r$（カリグラフ体）。$T_c$ は**セル $c$ の累積4次元アフィン変換**で別物。当初 $T_r$ と書いて衝突させた。

### JCDCG³ 発表スライド（`jcdcg3_slides.tex`，2026-08-19 仮作成）

**発表は 9/8（火）16:42--16:57**（プログラムで確認）。本編10枚＋**backup 5枚**，**4:3**（2026-09-03 に 16:9 から変更）。合計 12:30 なので**質疑は実質2分強**。

**構成の要点：Warning を単独の話題にせず，タイブレークの結果の系として置く。**

提出済みアブストラクトの Remark は「タイブレークが等変でないので全根を個別検証した。**全 BFS 木への拡張は未解決**」と自ら書いている。そこで

1. スライド7で**そのアブストラクトの一文を引用**して問題提起
2. 1,132,661本の検証で答えを出す（5/8/16/24胞体は規則に無関係）
3. スライド8で **⚠ Warning：一様には拡張しない**

の順にすると，追加情報ではなく**提出内容への回答**として聞こえる。逆順だと「一様の話のついでに正則も」となりアブストラクトから離れる。

| # | 内容 | 目安 |
|---|---|---|
| 1 | Title（**提出時のタイトルのまま**。プログラムとの整合） | 0:20 |
| 2 | 問題設定・D&H の3未解決ケース | 1:20 |
| 3 | Face-rotation BFS の3ステップ | 1:20 |
| 4 | **Step 2 hides a choice**：タイブレークの問いを立てる（2026-08-19 追加） | 0:50 |
| 5 | Validity test と `RegionMember` の罠 | 1:20 |
| 6 | 主定理＋全6胞体のネット（`260612AllFaceRotation.pdf`） | 1:40 |
| 7 | apple-peel との対比 | 1:30 |
| 8 | **And more**：タイブレークの問いに答えた＋**BFS 制限が空でないこと** | 1:50 |
| 9 | **⚠ Warning**：一様は規則に敏感 | 1:40 |
| 10 | Open problems ＋ GitHub URL（**24胞体は settled に訂正**） | 0:40 |
| B1 | backup：探索が到達できる木は 24胞体で 7.9% だけ | — |
| B2 | backup：3規則の比較表 | — |
| B3 | backup：runcinated 24-cell は等変な木なら valid（スライド4の伏線回収） | — |
| B4 | backup：セル次数も辺の角度欠損も両方向で効かない | — |
| B5 | backup：正多胞体6種の基本データ（論文 Table 1 ＋ 深さ・木の数） | — |

**合計 12:30**（2026-09-03 に 12:20 から。スライド8に Seamons の一文を追加した分）。枠は 16:42--16:57 なので質疑は約2分。**縮めたくなった場合の削りどころはスライド6（図を見せる時間 −10秒）とスライド2（−10秒）** で，スライド4は削らない。

**Warning は grand antiprism 一択**（`rule A → all 320 / rule B → no root`）。**64種の三分類（61/2/1）は出さない** — 数字自体が規則依存で留保が要り，切れ味を鈍らせるため。締めは "What regularity buys is not that some tree works, but that the choice does not matter."

#### スライド4（タイブレークの問い）の設計（2026-08-19 追加）

**発端**：査読でタイブレークが問われることが今の状況を生んだ以上，その問いを聴衆より先に自分で口に出すべき，というユーザー判断。

**スライド4では問いだけを立て，答えは一切出さない。** 締めは "We come back to this." で，答えはスライド8に温存する。順序を逆にしたり4で答えを漏らしたりすると，スライド8が「アブストラクトへの回答」ではなく単なる繰り返しになる。

構成は3点のみ：層は強制されるが親は強制されない／幾何的な規則に替えても救われない／だから規則は恣意的にならざるを得ない（我々のは cell index で，これも等変でない）。右段に $L_0$–$L_2$ の TikZ 図（$v$ に破線の候補親2本と「?」）。**`\usepackage{tikz}` を明示的に追加**（metropolis が内部で読むが依存しない）。

**主張は 5/8/16/24胞体に限定する。** 120/600胞体では成り立たない（上の訂正ボックス参照）。ちょうどスライド8で全木を尽くせる4種と一致するので話の筋も揃う。**例外はスライドに出さず「もし突っ込まれたら」としてノートに置く。**

#### スタイル確定（2026-08-19）

テーマを `default`+`dove` から **metropolis** に変更。

```latex
\usetheme[progressbar=frametitle,numbering=fraction,block=fill]{metropolis}
\usepackage[sfdefault]{FiraSans}   % pdflatex で Fira を使う（XeLaTeX 不要）
\usepackage{newtxsf}
\usepackage{appendixnumberbeamer}  % backup を「n/9」の分母から外す
```

- **採用理由は `block=fill`**：dove では alertblock の枠が出ず Warning スライドのインパクトが弱かった。metropolis の塗り潰しブロック＋`block title alerted` の赤指定で解決
- `navigation symbols` と `footline` の手動設定は metropolis が内蔵するため削除

**レイアウトの変更**（metropolis は余白・字送りが大きく，そのままでは Overfull が 3→6 箇所に増える）：

| スライド | 変更 |
|---|---|
| 1 Title | `[plain]`。短縮タイトル `\title[...]{...}` を追加（ノート頁で `\\` が詰まり "of theSix" になるため） |
| 5 Main result | 縦積み → **2段組**（左：定理ブロック＋773根，右：Fig 1）。図が大きく取れる |
| 6 Contrast | 表を `\small`，無題ブロックに `Consequence` の見出し（metropolis の塗り潰しでは空題が灰色バーになる） |
| 7 And more | 縦積み → **2段組**（左：アブストラクト引用＋説明，右：木の数の表）。結論行は全幅 |
| 8 Warning | 左段 `\small`，`talk_gap.png` を 0.72 幅に，締めの一文を `\small` で1行に |
| B2 | 表を `\small` |

**残る Overfull は title frame の 1件のみで，これは metropolis の構造的なもの**（title page を `minipage[b][\paperheight]` で組むため frame に収まりようがない）。出力は欠けていない。metropolis が pdflatex で出す "compile with XeLaTeX" 警告も FiraSans で代替済みのため無視してよい。両方ヘッダコメントに明記済み。

#### 8胞体の BFS ネットは常に Dalí の十字架（2026-08-19 確認）

スライド2の図を 5胞体から 8胞体（`face_rotation_net_8cell.png`）に差し替えた際に確認した事実。

**8胞体のセル隣接グラフは K_{2,2,2,2}**（対蹠セル以外はすべて隣接，全セルが次数6）。よって任意の根からの BFS 層は **1 + 6 + 1** で固定され，BFS 木は「対蹠セルを6個の隣接セルのどれにぶら下げるか」の6通りしかない（CLAUDE.md の「根あたり6本」と一致）。

展開後のセル重心を実測（辺長単位，根セルを原点）：

```
根 {0,0,0} ／ 6個の隣接セル {±1,0,0},{0,±1,0},{0,0,±1} ／ 対蹠セル {0,-2,0}
```

**＝ Dalí の十字架**（4個の柱＋中央セルの側面4個）。6本すべてで対蹠セルは「親と一直線上の距離2」に着地し，違いはどの軸が長腕になるかだけ（立方体の対称群で合同）。

つまり 8胞体では **48通り（8根 × 6木）すべてが Dalí の十字架と合同**。

#### 5胞体も全ネットが合同（2026-08-19 確認）

セル隣接グラフは **K₅**（全セルが互いに隣接，次数4）。層は **1 + 4** で層2が存在せず，各非根セルの親候補は根のみ → **根あたり BFS 木はちょうど1本**。したがって**ネットは全部で5個しかない**。

形状：**中央の正四面体の4面すべてに正四面体を貼った形**（相異なる頂点は8個）。外側4セルの重心は根の重心から 1/√3 ≈ 0.5774，互いに 2√2/3 ≈ 0.9428 で正四面体をなす。

**5個すべてが合同**：頂点間距離の多重集合が 1.3×10⁻¹⁵ で一致し，Kabsch で明示的に重ねると残差 ~5×10⁻¹⁶。**det = +1 の真の回転で重なる**（鏡映不要。形自体が achiral なので det = −1 の解も同時に存在するが，混同しないこと）。

#### 合同性はここで止まる

| | セル | 次数 | 根1からの層 | 根あたり BFS 木 |
|---|---:|---:|---|---:|
| 5胞体 | 5 | 4 | 1 + 4 | **1** |
| 8胞体 | 8 | 6 | 1 + 6 + 1 | **6**（全て合同） |
| 16胞体 | 16 | 4 | 1 + 4 + 6 + 4 + 1 | 20,736 |
| 24胞体 | 24 | 8 | 1 + 8 + 14 + 1 | 32,768 |

5胞体・8胞体で全ネットが合同なのは**層が浅く選択の余地が無い／あっても対称で吸収される**ため。16胞体以降は層が深く親候補が実質的に分岐するので，同じ議論は使えない（全木 valid であることとネットが合同であることは別問題）。

検証スクリプトはスクラッチのみで未保存。`src/bfsnet/face_rotation_net_all4D_v2.m` を `Print["Face-rotation BFS net check v2` の直前で切って `ToExpression` し，`faceUpForRootP` / `bfsUnfoldP` / `applyAff` で重心・頂点を出すだけで再現できる。

#### スピーカーノート（2026-08-19 追加）

本編全10枚に `\note{}` を追加。**目標時間（個別＋累計）を各ノートの先頭に太字で置く**（合計 **12:20**）。プリアンブルのコメントを外すだけで2種のビルドが得られる（両方コンパイル確認済み）：

```latex
\usepackage{pgfpages}
\setbeameroption{show notes on second screen=right} % 発表用．907×255pt の横長1枚に
\setbeameroption{show only notes}                   % ノートのみ．印刷用
```

- 通常ビルドではノートは一切出力されない（12ページのまま）
- 2画面版は Skim / Adobe Reader の「プレゼンモード＋2画面」で使う。ページ番号が `8/100` と誤表示されるが**ノート版のみの症状**で本番 PDF には影響しない
- backup B1/B2 にはノートなし（質疑応答用のため）

#### 図の背景を統一（2026-08-19，解決）

当初「Fig 1 は黒背景グレースケール，5胞体と grand antiprism は白背景カラー」という混在があった。2段階で解消：

1. スライド2の図を 5胞体 → **8胞体**（`face_rotation_net_8cell.png`）に差し替え
2. `src/figures/make_talk_figures.m` の `renderNet` に **`Background -> GrayLevel[0.13]`** を追加し `talk_gap.png` を再生成
3. `src/bfsnet/face_rotation_net_viz_all.m` の背景を `Black` → **`bgCol = GrayLevel[0.13]`**（冒頭で定義）に変更し 8胞体を再生成

**どちらもセルの配色は変更していない**（gap は BFS 深さのグラデーション，8胞体は GrayLevel の濃淡）。

実測した背景値（グレースケール 0–255）：

| 図 | スライド | 背景 |
|---|:---:|---:|
| `face_rotation_net_8cell.png` | 2 | **33** |
| `260612AllFaceRotation.pdf` の各パネル | 6 | **31** |
| `talk_gap.png` | 9 | **33** |

31 と 33 の差は不可視。**スライドの3枚は揃った。**

> **`src/bfsnet/face_rotation_net_viz_all.m` の PlotLabel について**：このスクリプトの `PlotLabel` は元々 `Black` で背景も `Black` だったため**一度も見えたことがない**（実測：上端60行が全て輝度0）。背景だけ変えると黒い幽霊文字が浮くので，**ラベル色も `bgCol` に変更**して不可視のまま維持した。削除しなかったのは，ラベルが上端の帯を確保しており，消すと 600×600 のキャンバス内でネットが再センタリング・拡大されて**構図が変わってしまう**ため。スライド側には独自のキャプションがある。

> **未再生成**：`face_rotation_net_{5,16,24,120,600}cell.png` は**純黒（0）のまま**。スライドでは使っていないので実害はないが，スクリプトの既定は 0.13 になったので次に流せば揃う。600胞体は hull 構築が重い。

> **注意**：`260612AllFaceRotation.pdf` の**ページ角は白**（255）。CLAUDE.md が言う GrayLevel[0.13] は**パネル内部**の値で，角をサンプルすると白が返る。

#### 24胞体が all-net でないことの反映（2026-09-03）

[[project-24cell-not-allnet]]（`memory/project_24cell_not_allnet.md`）の結果をスライドに反映。
**スライド10の2つ目の項目が事実として誤りになっていた**ため訂正は必須だった。

| スライド | 変更 |
|---|---|
| 10 Open problems | 「still open for the 24-, 120- and 600-cell」→ **「settled for the 24-cell (Seamons 2026), still open for the 120- and 600-cell」**。ノートは「最初と最後だけ読む」なので**発話時間は増えない** |
| 8 And more | 最終行に一文追加（+10秒）：*"And BFS is not a vacuous restriction: the 24-cell is **not** all-net (Seamons 2026) --- some spanning tree does overlap."* |

**単なる訂正ではなく講演が強くなる**：「なぜ BFS 木に限るのか，コードがたまたまそうだっただけでは」という当然の反論に答えられるようになる。**24胞体では一般の全域木に重なるものがあるのに BFS 木は 786,432本すべて valid** ＝ BFS はちょうど境界線であって恣意的な制限ではない。スライド9の締めと同じ構図を正則多胞体の内部でもう一度見せることになる。

スライド8のノートに，この一文が「why BFS?」への回答であること（弁解ではなく理由として言うこと）と，**出典を突かれた場合の答え**（2026年8月の Zenodo preprint，doi:10.5281/zenodo.22004769，自分で厳密有理数で再現済み）を追記。

> **著者名の綴りに注意：Seamons**（Seamos ではない）。

#### 4:3 への変更（2026-09-03，完了）

**会場のプロジェクタが 4:3**（ユーザーが 2024 年の JCDCG³ で確認）。`aspectratio=169` → **`43`**。

4:3 は 16:9 より**横が狭い**（128mm 対 160mm）だけで縦は逆に高い。段が細くなると同じ本文でも行数が増えるので，**2段組の frame が軒並み縦にあふれた**（frame 2 / 4 / 8 / 9 の4枚）。とくに frame 4 は 49.5pt 超過。

| frame | 対処 |
|---|---|
| 2 Unfolding | 左段を `\small`，`\medskip`→`\smallskip` |
| 4 Step 2 hides a choice | 段を **0.57/0.39 → 0.68/0.29**，tikz `scale` 0.95 → **0.62**。`\small` は維持 |
| 8 And more | 左段 `\small`，引用を `\raggedright`（underfull 解消），表は `\small` のまま `\tabcolsep` を **2pt**（表が 5.5pt はみ出していた），下段の総括文を `\small` |
| 9 Warning | 段を 0.55/0.42 → **0.58/0.39**，図 0.72 → **0.86**（段が細くなったぶん），`\smallskip` を1つ削減 |

**frame 4 は font size を落とさずに済んだ**。`\footnotesize` にすれば通るが，段幅を 0.68 に広げるだけで `\small` のまま収まる（tikz の `scale` は座標だけを縮め，ノードの文字サイズは変わらないので図は密になるが破綻しない）。

> **落とし穴**：frame 9 で `\smallskip` を削るとき**空行ごと**消すと段落の切れ目まで消え，「Same polytope…」と「Among the 64…」が**1段落に融合する**。空行は残して `\smallskip` だけ削ること。

**残る警告はタイトル枚の overfull vbox 1件のみ**（metropolis が title page を `\paperheight` の minipage で組む構造上のもの，従来どおり無害）。**12ページ**維持。

段幅・`\small`・tikz scale は**すべて 4:3 前提で調整済み**なので，16:9 に戻すとコンパイルは通るが本文段が不必要に細いままになる。この経緯は `jcdcg3_slides.tex` のヘッダコメントにも記載。

#### 2026-09-06（発表2日前）：動画・印刷用・ノートの上限

**引用の裏取り**：スライド8が引く *"Whether validity extends to all BFS tie-breakings ... remains open."*
は**採択後の camera-ready `260702JCDCGGG_FaceRotationNet.pdf` にある**（初回提出の 0611 版には無い）。
省略した後半 "and further to all spanning trees" もスライド8末尾の Seamons とスライド10で答えており，
**引用文の前半と後半を講演が両方埋める**構図になっている。

**修正3件**：(1) B2 の「61種のうち13種，その中に2640セルの4種」→ 2640セルの4種のうち**2種は MIXED で
61に入っていない**（population 違い）。同じ誤りが本ファイルと論文§4にもあった。(2) スライド9
`others swing` → `another swings`（1.8→99.8% は runcitruncated 600-cell 1種）。(3) 質疑時間を2分に訂正。

##### 埋め込み動画（3本，pympress で再生）

`\usepackage{multimedia}` の `\movie` が書く PDF Movie 注釈を pympress が Poppler 経由で拾い，
GStreamer で再生する（`avdec_h264` と Apple の `vtdec` を実機確認）。**2画面版でも注釈の矩形が
正しく半分に縮む**ことを Poppler で実測確認済み。

| スライド | ファイル | 内容 |
|---|---|---|
| 2 | `talk_8cell_spin.mp4` | 完成したネットのターンテーブル。10秒で1回転 |
| 3 | `talk_8cell_unfold.mp4` | テッセラクト → Dalí の十字架 → 再び折り畳む。**層ごとに開く** |
| 6 | `talk_mainresult_spin.mp4` | Figure 1 の**下段2枚だけが同時に回る**。8秒で1回転 |

**poster に静止画を置くのが肝**：メディア非対応のビューアでは注釈が無効になり poster が出るので，
**素の4:3 PDF を会場PCに渡しても一切劣化しない**。**クリックで再生開始**（Poppler が autoplay
フラグを公開しないため pympress に自動再生はない）。**mp4 は PDF と同じディレクトリに置く**。
生成は `./build_talk_video.sh`（`src/figures/make_talk_video_{8cell,unfold,mainresult}.wls` ＋ ffmpeg，約5分）。

**展開アニメ（スライド3）**：$T_r(t)=\mathrm{id}$，$T_c(t)=\mathrm{rot}(t\theta_e)\circ T_p(t)$。
蝶番平面は親の**現在の**フレームで取り直し，$\theta$ は子を**まだ回していない**位置に対して測るので
$t$ に依らない。`coords4` は**層ごとの展開率のベクトル**を取り，6個 → **0.6秒の静止** → 最後の1個
という順に開く。素直に $w$ を落とすと $t=0$ でただの立方体になるので，$w$ 軸上の視点による
**4次元透視投影**を使い，重みを $1-t$ で消して $t=1$ でちょうど plain drop-$w$ ＝本物のネットにする。

> **描画の落とし穴**
> - 透視を最後まで効かせると，$w=-2.5$ まで降りるセルが縮み**広がる前に一度小さくなる**
> - **`SphericalRegion -> True` は展開アニメでは外す**（外接球に合わせるので $\sqrt3$ 損する）。
>   カメラが動く turntable では必要
> - `PlotRange` は**クリップする**ので幾何より狭められない。しかも直方体の投影は六角形なので
>   外接矩形との差で**3分の1が空く**。→ 大きめに描いて **`ffmpeg` の `cropdetect` で全フレームの
>   描画範囲の和を測ってクロップ**（`build_talk_video.sh` の `tight_crop`）。数値をハードコード
>   しないのでカメラを変えても自動追従する
> - スライド6は `260612AllFaceRotation.pdf` が**Illustrator の組版**なので再現不能。元の図を
>   100 dpi で起こし，**下段2枚だけを同じ矩形に overlay** する（562×560，y=302，x=12 と 583）。
>   元パネルが `face_rotation_net_viz_all.m` 由来と判明したので，同じ `ViewPoint {2.4,-2.0,1.8}`
>   から回し始めればフレーム0が印刷図と一致する

##### ノートが切れるのは LaTeX の問題（pympress ではない）

beamer のノート頁は**改頁しない単一の vbox** なので，はみ出した分は頁の外へ押し出され **PDF に
入らない**。`pdftotext` で確認したところ，スライド2・3・8のノートは末尾数行が PDF に存在しなかった。

> **ログは当てにならない**：TeX は overflow が大きいときしか `Overfull \vbox` を出さない。
> 3件切れていたのに警告は2件。しかも `build_slides.sh` は素のビルドのログしか見ていなかった
> （素のビルドではノートを組まないので永久に気づけない）。

**`check_notes.py`** を追加（`build_slides.sh` から自動実行）。各 `\note` の末尾5語を，**ノート領域
だけを `pdftotext -x/-y/-W/-H` で切り出した文字列**と照合する。領域を限るのは必須で，スライド本体と
同じ語で終わるノート（スライド4の "we come back to this"）は頁全体を見ると**切れていても通る**。
2026-09-06 に全10枚のノートを短く平易な英語に書き直した（目安は組み上がり12行）。

##### 3種のビルドと印刷用

ノートの切り替えはコメントの付け替えをやめ **`\ifdefined\NOTES` / `\ifdefined\NOTESPRINT`** に変更。
`./build_slides.sh` が3種を焼き，警告数・ノートの完全性・持ち物を表示する。

- `jcdcg3_slides.pdf` 素の4:3（会場PC用。**必ず持参**）
- `jcdcg3_slides_notes.pdf` 横長2画面版 → `pympress -t 12:30 -N right ...`
- `jcdcg3_slides_print.pdf` **印刷用**。`show only notes` ＋自前の `note page` テンプレートで
  **左にスライド・右にスクリプト**，`\note` を持つ frame だけなので10ページ

> **印刷用テンプレートの落とし穴**：幅は **`\paperwidth` ではなく `\textwidth` の割合**で指定する
> （ノート頁もスライドの余白を保つ）。`\insertslideintonotes{f}` の返す箱は **f×paperwidth より広い**
> （係数およそ 388pt）ので，左の minipage をそれより広く取らないと overfull hbox が10ページ分出る。
> 現行は scale 0.40／minipage 0.48・0.49。本文は `\footnotesize`（`\small` だとスライド4と8が切れる）。

**pympress 1.8.6**（homebrew，導入済み）。既定の `[notes position] horizontal = right` が
`show notes on second screen=right` と一致する。キー：`s` 画面入れ替え，`b` 暗転，`n` ノートモード，
`p` タイマー一時停止，`r` リセット，`f` 全画面。

### 「良いネット」の三層定義と対称性の上限（2026-08-19）

**発端**：ユーザーの問い「良い net とは何か？　私は対称性が良いと考える」。Wikipedia 英語版の uniform 4-polytope のネットが validity を考慮していないのに対し，3D の展開図では validity が必須であることは自明 → **4D でも validity は必須**というのがユーザーの立場。

#### 定義（ユーザー承認済み，2026-08-19）

| 層 | 基準 | 位置づけ |
|:--:|------|---------|
| 1 | **validity**（内部が重ならない） | **必須条件** |
| 2 | **対称性 \|H\|** の最大化 | 主基準 |
| 3 | **コンパクトさ** | 副基準 |

#### 層2：対称性は木の選び方の言葉に翻訳できる

face-up 後，$\mathrm{Stab}(r)$ は $w$ 軸を固定するので $xyz$ 超平面上の $O(3)$ 部分群として作用する（＝根セルを3次元多面体として見た対称群）。展開の各蝶番回転は幾何学的に定義されるので $T_{h(c)} = h \circ T_c \circ h^{-1}$ が成り立ち：

> **$H \le \mathrm{Stab}(r)$ が BFS 木 $T$ を根付き木として保つなら，ネットは $H$ 対称。**
>
> **$H$-不変な BFS 木が存在 $\iff$ すべてのセル $v$ について親候補 $C(v)$ が $\mathrm{Stab}_H(v)$ の固定点を含む。**

実行可能な部分群の族は**下に閉じている**（$K \le H$ が feasible なら $H$ の固定点はそのまま $K$ の固定点）ので，feasible な巡回部分群から生成元を足して登る探索で**最大値が厳密に求まる**。

**これはタイブレーク問題への回答でもある**：「候補を等変に区別する」のではなく「**区別しないまま等変に配る**」。恣意性は各軌道につき1つの代表選択に縮む。

#### 対称性の上限（全6種，確定）

`src/symmetry/symmetry_ceiling.wls`（根 r=1，Aut(セル隣接グラフ) を全列挙）

| | 根セル | \|Stab(r)\| | 現行(queue) | **上限 MAX \|H\|** | 最大部分群の個数 |
|---|---|---:|---:|---:|---:|
| 5胞体 | 正四面体 | 24 | **24** | **24** | 1 |
| 8胞体 | 立方体 | 48 | **8** | **8** | 3 |
| 16胞体 | 正四面体 | 24 | 1 | **1** | — |
| 24胞体 | 正八面体 | 48 | **6** | **6** | 4 |
| 120胞体 | 正十二面体 | 120 | **1** | **10** | 6 |
| 600胞体 | 正四面体 | 24 | 1 | **1** | — |

**現行実装が最適でないのは 120胞体だけ**（1 → 10 の改善余地）。他の5種は既に上限を達成している。

#### 経験則：上限は「蝶番1枚の対称群」

立方体・正八面体・正十二面体セルでは

$$\text{MAX}\,|H| \;=\; \frac{|\mathrm{Stab}(r)|}{(\text{根セルの2-面の数})} \qquad 48/6=8,\quad 48/8=6,\quad 120/12=10$$

が**3種とも厳密に一致**し，最大部分群の個数も（面の数）/2 = 3, 4, 6 と一致する（対蹠面の安定化群は同一）。解釈：**ネットはどこかで根セルの2-面を1枚区別せざるを得ず，残るのはその面を保つ対称性だけ**。

8胞体では証明できる：対蹠セルは $\mathrm{Stab}(r)$ 全体に固定されるのに親候補6個が1軌道をなす → $H$ は6個のうち1個を固定 → $|H| \le 48/6 = 8$，かつ達成。**24/120胞体の一般証明は未着手。**

**正四面体セルは例外で全崩壊**（16胞体・600胞体で MAX = 1）。16胞体（セル隣接グラフ $= Q_4$）では単位元以外のどの元も単独で不可能：$v = 1100$ は互換 $(12)$ に固定されるのに親候補 $1000, 0100$ が入れ替わる。3-巡回・4-巡回・二重互換でも同型の障害。**つまり 20,736 本の BFS 木は1本残らず非対称。** 600胞体も同様に MAX = 1。「対称性が良さである」という基準が**全域的に消える多胞体が存在する**。

#### 層3：外接球半径は使えない（重要）

ユーザー提案は「外接球半径」だったが**飽和する**。無作為な BFS 木8本で測った範囲：

| | R_minball | R_gyration | 凸包体積 |
|---|---|---|---|
| 8胞体 | 4.242641（一定） | 2.783882（一定） | 144（一定，全木合同なので当然） |
| 16胞体 | 2.2233–2.3246 | 1.3180–1.3467 | 13.68–14.46 |
| 24胞体 | 5.033223（**一定**） | 3.0225–3.0338 | 253.5–257.2 |
| 120胞体 | 8.531294（**一定**） | 4.7473–4.7642 | 1445–1450 |

**24胞体と120胞体では最小外接球の半径が木に依らず完全に一定**（15桁一致。中心だけが動く）。理由は未解明。よって**副基準は回転半径（radius of gyration）を採る**。全ネットの体積は同じ（合同なセル $n$ 個）なので公平な比較になる。

> **留保**：回転半径と凸包体積は**順位が一致しない**。120胞体の $H$-不変ネットは無作為な木より回転半径が小さい（4.7224 < 4.7473）が凸包体積は大きい（1463 > 1450）。「質量の広がり」と「必要な外接空間」は別物。どちらを採るかは設計判断。

#### 120胞体：扇形構想は最適だった（`src/symmetry/sector120_nets.wls`）

ユーザーが 2026-08-17 に発案した「根の正十二面体の五角形面の軸で5回対称に分割する」案は，**計算上ちょうど最適**だった。

- $H = D_5$（位数10，五角形面の安定化群），そのような軸が **6本**（12面の対蹠対）
- $H$-軌道は **31個**（サイズ1が10個，5が20個，10が1個）
- $H$-不変な BFS 木は **2048本**（自由度のある軌道が9個 + 4択が1個）
- **2048本すべて valid**（重なり0）→ **対称性と validity は競合しない**
- 回転半径 4.722389 〜 4.781145（相異なる値 563個），**最小は2本ちょうど**が達成（鏡像対と思われる）
- 最良木：index 1824，R_gyr = 4.722389，凸包体積 1463.28

データは `data/sector120_setup.mx` / `data/sector120_best.mx`。

#### 16胞体・600胞体に対称ネットが無い理由（2026-08-19，決着）

**ユーザーの問い**：「対称に伸ばしても invalid になるのか？　それとも木の構造そのものが対称にならないのか？」

**答えは後者。しかも予想より強い。validity は一切関与しない**（計算はすべて組合せ的で，幾何を見ていない）。

**補題（証明済み）**：根付き全域木 $T$ が $H$ 不変なら，$v \in \mathrm{Fix}(H)$ の親も $\mathrm{Fix}(H)$ に属する（親は $\mathrm{Stab}_H(v) = H$ の固定点でなければならないため）。したがって

> **$\mathrm{Fix}(H)$ はセル隣接グラフの中で，根を含む連結部分グラフを誘導しなければならない。**

**BFS を外しても変わらない**（`src/symmetry/symmetry_spanning_trees.wls`）。親を「1層上」に限らず任意の隣接セルに許した一般の全域木でも，上限は 24 / 8 / **1** / 6 / 10 / **1** で BFS と完全に同一。**障害は最短路条件ではない。**

`src/symmetry/symmetry_fixed_column.wls` による $\mathrm{Fix}(h)$ の層分布と連結成分数（各元位数の代表1個）：

| | 位数 | \|Fix\| | 層ごとの分布 | 成分数 |
|---|---:|---:|---|---:|
| 16胞体（層 1,4,6,4,1） | 2 | 8 | {1,2,2,2,1} | **2** |
| | 3 | 4 | {1,1,**0**,1,1} | 2 |
| | 4 | 2 | {1,0,0,0,1} | 2 |
| 600胞体（層16段） | 2 | 60 | 途中に 0 あり | **12** |
| | 3 | 12 | 途中に 0 多数 | 6 |
| | 4 | 2 | 根と対蹠のみ | 2 |
| **120胞体** | **2** | **30** | {1,4,8,8,8,1} | **1** ✓ |
| | **5** | **10** | **{1,2,2,2,2,1}** | **1** ✓ |
| | 3 / 6 / 10 | 6 / 2 / 2 | 途中に 0 | 6 / 2 / 2 |

**16胞体・600胞体では非単位元 23個すべてで $\mathrm{Fix}(h)$ が非連結。** よって $|H| = 1$ が全域木レベルで確定。

**幾何学的な理由（回転の場合）**：固定セルの連なりは，軸が**面から面へ**セルを貫くときにだけ続く。ところが**正四面体の対称軸は面と面を結ばない** — 3回軸は面↔対頂点，2回軸は辺↔対辺。だから正四面体セルの多胞体（16胞体・600胞体）では軸が1歩でセルの外に出てしまい，柱が切れる。16胞体の3回軸で**層2の固定セルが0個**になるのがその直接の現れ。

対して立方体・正八面体・正十二面体は**心対称で軸が対面を貫く**ので柱が続く。120胞体の5回軸が固定する10セル {1,2,2,2,2,1} が，まさにユーザーの言う「扇形の軸」である。**これが「上限＝根セルの2-面1枚の安定化群」の正体**でもある。

> **留保**：鏡映（位数2）は軸の議論では説明できない。16胞体では鏡映の固定セルが全層に存在するのに2成分に割れる（$Q_4$ で $\mathrm{Fix}((12)) = \{(a,a,c,d)\}$ が $a=0$ と $a=1$ の2ブロックに分かれ，両者はグラフ距離2で辺が無い）。連結性の判定は計算で確認済みだが，統一的な証明は無い。

> **なお，「対称だが invalid」の例はまだ1つも無い。** 対称な木が存在する4種（5/8/24/120胞体）では，調べた限り全て valid（120胞体の $H$-不変木 2048本すべて）。探すなら一様多胞体の MIXED 2種・ALL INVALID 1種。

#### runcinated 24-cell は ALL INVALID ではなかった（2026-08-19，重要）

**`x3o4o3x`（runcinated 24-cell）は「凸一様64種で唯一の ALL INVALID」と記録していたが，これは3つの index 系規則についての事実であって多胞体の性質ではない。対称な木を課すと valid なネットが出る。**

セル隣接グラフの $|\mathrm{Aut}| = 2304 = 1152 \times 2$（24胞体の自己双対性による拡大 F₄ 群）。セル軌道は2つ（八面体48・三角柱192）。どちらの軌道の根でも **MAX \|H\| = 3**（$C_3$）。

`src/symmetry/uniform_symmetry_nets.wls`，各600本サンプル：

| 根 | 軌道 | \|Stab(r)\| | MAX \|H\| | **$H$-不変な木** | **無制約ランダム木（対照）** |
|---|---:|---:|---:|---|---|
| 1 | 48（八面体） | 48 | 3 | **VALID 70/600**（中央値6） | **VALID 0/600**（中央値7） |
| 25 | 192（三角柱） | 12 | 3 | **VALID 83/600**（中央値6） | **VALID 0/600**（中央値9） |

**対照実験が本質。** 無制約のランダム木では 1200本中0本が valid なのに，3回対称を課すと 12〜14% が valid になる。**対称性そのものが効いている**（規則が index でないことの効果ではない）。$H$-不変な木では重なり数が必ず3の倍数になる（重なりが $C_3$ 軌道をなすため）。

> これは「validity を必須とした上で対称性を最大化する」という三層の設計が，**単なる美的な選好ではなく valid なネットを見つける手段**でもあることを示す最初の実例。

#### MIXED 2種（H₄，2640セル）：対称性は常に効くが万能ではない

セル軌道は両種とも4つ（120 / 600 / 720 / 1200）。軌道ごとに根を1つ取り，$H$-不変な木と無制約ランダム木を各25本ずつ。

**`x5o3x3x`（runcitruncated 600-cell，queue 規則で 47/2640）**

| 根 | 軌道 | \|Stab(r)\| | MAX\|H\| | $H$-不変（valid / 最小重なり / 中央値） | ランダム（同） |
|---|---:|---:|---:|---|---|
| 2521 | 120 | 120 | 5 | 0/25 ・ 15 ・ 35 | 0/25 ・ 23 ・ 43 |
| 1 | 600 | 24 | **1** | 0/25 ・ 35 ・ 45 | 0/25 ・ 35 ・ 45 |
| 1801 | 720 | 20 | 5 | 0/25 ・ **5** ・ 20 | 0/25 ・ 16 ・ 24 |
| 601 | 1200 | 12 | 2 | 0/25 ・ 12 ・ 26 | 0/25 ・ 21 ・ 32 |

**`x5x3o3x`（runcitruncated 120-cell，queue 規則で 2570/2640）**

| 根 | 軌道 | \|Stab(r)\| | MAX\|H\| | $H$-不変 | ランダム |
|---|---:|---:|---:|---|---|
| 1 | 600 | 24 | **1** | 0/25 ・ 4 ・ 12 | 0/25 ・ 7 ・ 12 |
| 1801 | 720 | 20 | 5 | **5/25** ・ 0 ・ 5 | 1/25 ・ 0 ・ 6 |
| 601 | 1200 | 12 | 3 | **7/25** ・ 0 ・ 3 | 1/25 ・ 0 ・ 5 |
| 2521 | 120 | 120 | 5 | **6/25** ・ 0 ・ 10 | 0/25 ・ 7 ・ 12 |

**読み方**

- **対称性が使える全ケースで重なりの中央値が下がる。悪化した例はゼロ。** valid 率も上がる（x5x3o3x で 5〜7/25 対 0〜1/25）
- **ただし x5o3x3x は 25本では1本も valid にならない**（最良は根1801の重なり5個＝$C_5$ 軌道1個ぶん，valid 直前）。対称性は効くが決定打にならない
- **MAX\|H\| = 1 の根（両種の root 1）では2行が同じ分布になる** — 制約が消えるので当然であり，**パイプラインの内部整合性チェックになっている**（x5o3x3x root 1 は min/median が 35/45 で完全一致）
- 一様多胞体では**根によって MAX\|H\| が 1〜5 と大きく変わる**。正多胞体（セル推移的なので根に依らない）と異なり，**根選びが対称性の上限を決める**

> **サンプリングの意味**：「$H$-不変な木上の一様分布」対「全 BFS 木上の一様分布」の比較。既存の `47/2640` 等は**根ごと・規則固定**の数字なので直接は比較できない。同一根・同一本数での対称 vs 無制約の対比が本実験の主張。

### 辺まわりの角度欠損 — 空間充填仮説とその反証（2026-08-19）

**ユーザーの発案**：「胞を作る多面体が空間充填できるかは指標になりうる。いや，充填できなくても面を重ね続けたときに変な曲がり方をしなければいいのかな」。

**この2つの言い方は同じ量を指す。** `src/symmetry/edge_angular_defect.wls`：

$$\mathrm{defect}(e) \;=\; 2\pi - \sum_{c \ni e} (\text{セル } c \text{ の } e \text{ における二面角})$$

- **defect = 0 ⟺ セルが辺のまわりにちょうど収まる**＝局所的に R³ の面接触タイル張り。凸性から常に defect > 0 なので**凸4-多胞体は決してタイル張りにならない**が，いくらでも近づける
- **defect は展開のホロノミー**でもある：辺のまわりを1周すると，合成された面回転は defect ぶんの回転になって戻る。つまりネットは**辺に沿って錐角 $2\pi - \mathrm{defect}$ の錐特異点を持つ錐構造の展開写像**。「変な曲がり方」の正体
- **多胞体が Euclid 蜂の巣からどれだけ離れているかの尺度**：{4,3,3} = 90°，{4,3,4}（立方体蜂の巣）= 0，{4,3,5} < 0（双曲的）

#### 正多胞体6種：難しさの順序と完全に一致（辺は全て同値なので値は1つ）

| | 辺あたりのセル | defect |
|---|:--:|---:|
| 5胞体 | 3 | **148.41°** |
| 8胞体 | 3 | 90.00° |
| 16胞体 | 4 | 77.88° |
| 24胞体 | 3 | 31.59° |
| 120胞体 | 3 | 10.30° |
| 600胞体 | 5 | **7.36°** |

D&H が全スパニング木で解決した3種（5/8/16胞体）が上位3つ，apple-peel が失敗し全数列挙も不可能な 120/600胞体が下位2つ。**大きい defect ＝ 疎で余裕がある ＝ 展開しやすい。**

#### しかし一様多胞体では反証される（両方向で）

生成済み50種の最小 defect を測った結果：

| symbol | セル | min defect | 分類 |
|---|---:|---:|---|
| **x5x3o3x** | 2640 | **0.736°** | **MIXED 2570/2640** |
| **x5x3x3x** | 2640 | **0.736°** | **ALL VALID** ← 同値なのに valid |
| x5x3x3o / x5o3x3o | 1920 | 1.434° | ALL VALID |
| **x5o3x3x** | 2640 | **1.434°** | **MIXED 47/2640** |
| o5x3o3x / o5x3x3x / x5o3o3x | 1440–2640 | 1.471° | ALL VALID |
| … | | | |
| s3s4o3o | 144 | 10.224° | ALL VALID |
| x5o3o3o（120胞体） | 120 | 10.305° | ALL VALID |
| **x3o4o3x** | 240 | **10.529°** | **ALL INVALID** ← 高いのに invalid |

- **低い defect ⇏ invalid**：`x5x3x3x` は MIXED の `x5x3o3x` と**最小 defect が完全に同値（0.736°）**なのに全根 valid
- **invalid ⇏ 低い defect**：`x3o4o3x` は 10.529° で，多くの ALL VALID 種より高い

**よって角度欠損は「危険因子」ではあっても求める不変量ではない。** セル次数のときと同じく両方向に効かない（CLAUDE.md「セル次数は両方向に効かない」参照）。

> **「セルの多面体が空間充填するか」自体も指標にならない**：立方体（充填する）を持つ8胞体は易しいが，正四面体（充填しない）の5胞体はもっと易しい。逆に `x3o4o3x` は三角柱（充填する）を192個持つのに最難。効くのは**個々のセルの性質ではなく，辺のまわりの配置**である。

#### 積算量（defect × 深さ）も駄目（2026-08-19，`src/symmetry/curvature_battery.wls`）

「defect は局所量なので，木の深さで積算すれば巨大種と小さい種を同じ土俵に乗せられるのでは」という案を検定した。全50種について辺 defect の全リスト，総曲率 $\sum_e \mathrm{defect}(e)\cdot|e|$，セル隣接グラフの radius / diameter を計算し（`data/curvature_stats.mx`），**13個の候補統計量**で3つの問題種を分離できるか調べた。

| 統計量 | 3種の順位（昇順，全50中） |
|---|---|
| minD（生の最小欠損） | {1, **24**, 4} |
| totCurv / nC^(2/3) | **{6, 8, 10}** ← 最良 |
| (totCurv/nC) × diam | {5, 11, 13} |
| totCurv / nC | {4, 7, 15} |
| minD × diam | {4, 5, 29} |
| meanD × diam | {3, 8, 25} |
| diam / radius 単独 | {39,46,48} / {41,46,47}（逆向き） |

**どれも分離しない。** 最良の `totCurv/nC^(2/3)` でも3種は順位 6, 8, 10 で，**最下位2つは grand antiprism と snub 24-cell（どちらも ALL VALID）**，3位が600胞体（ALL VALID）。

さらに悪いことに，**積算すると生の defect が持っていた唯一の綺麗な信号が壊れる**。`totCurv/nC^(2/3)` での正多胞体6種の順位は 5胞体34 / 8胞体43 / 16胞体17 / 24胞体24 / 120胞体35 / **600胞体3** で，難しさの順序が完全に散らばる。

> **結論：曲率系のスカラー量はこの問題の不変量ではない。** セル次数・辺の角度欠損・曲率密度・曲率×深さ — **すべて両方向で反例が出る**。validity は単純な計量的密度の関数ではなく，展開写像が実際にどう折り返すかという**大域的・組合せ的**な性質だと考えるべき。

#### 未検証の仮説：「充填するか」はスカラーではなく二値として効くかもしれない（2026-08-19 記録のみ）

defect への翻訳で捨てた情報がある。**「セルが空間充填するか」と「defect が 0 か」は別物**である。

実例：8胞体のセルは立方体で空間充填するが，defect は 90°（0 ではない）。それでも展開後のセル重心は（辺長を単位として）

```
{0,0,0}, {±1,0,0}, {0,±1,0}, {0,0,±1}, {0,-2,0}
```

と**すべて整数格子点**に落ちる。8胞体の展開写像は **Z³ 上のウォーク**そのもの。

**仮説**：充填性はスカラーではなく**質的な二分**として効くのではないか。

| セルの性質 | 展開写像の像 | 衝突判定の性格 |
|---|---|---|
| 空間充填する（立方体・三角柱・切頂八面体…） | 離散集合（タイル位置） | **組合せ的**。「同じタイルに2つ載るか」だけ。ニアミスが原理的に存在しない |
| 充填しない（正四面体・正八面体・正十二面体…） | 無理数的な角度で降り積もる | **計量的**で微妙。「ぎりぎり当たる／外れる」が生じる |

120胞体で apple-peel が全滅するのに BFS では通るという繊細さは，正十二面体が充填しないことの現れかもしれない。

> **重要な留保**：defect が完全に 0 の理想的な格子でも validity は保証されない。そのとき問題は「木の全セルが相異なるタイルに載るか」という純粋な単射性の問題として**そのまま残る**。つまり充填性は問題を*離散化*するだけで*解決*はしない。
>
> **既に雲行きは怪しい**：`x3o4o3x` は三角柱（充填する）を192個持ちながら最難。全セルが充填多面体である種を集めても valid 率と対応しない可能性は高い。

**検定するなら**：64種を「全セルが充填多面体／一部／どれも充填しない」で三分し valid 率と突き合わせる。充填する凸多面体は有限リストで判定可能，計算は軽い。**未着手（ユーザー方針 2026-08-19：可能性として記録するにとどめる）。**

#### 凸一様64種の対称性上限（2026-08-19，`src/symmetry/symmetry_ceiling_uniform.wls`）

論文2の骨格を固めるため，**全64種 × 全セル軌道**で上限を計算。セル軌道ごとに根を1つ取れば十分（対称群が軌道に推移的なので同じ軌道の根は等価）。結果は `data/ceiling64.mx`，所要約25分。

**64種・170セル軌道。**

| 上限 \|H\| | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 | 24 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 軌道数 | **28** | 10 | 21 | 19 | 17 | **42** | 17 | 14 | 1 | 1 |

- **170根のうち28根（16%）が対称ネットを全く持たない**
- **どの根からも持たない種は3つ**：16胞体，600胞体，そして**runcinated 5-cell `x3o3o3x`（新発見）**
- **全安定化群を達成する種は2つ**：5胞体（24），そして**正四面体角柱（12，側面セル4個の軌道で）**（新発見）
- **上限が根に依らないのは 64種中18種だけ。** 最大の振れ幅は grand antiprism の **10 : 1**，次いで `x4x3o3o` と `o4x3o3o` の 8 : 1。正多胞体はセル推移的なので自明に一定

**runcinated 5-cell も同じ機構**（`Fix(h)` 非連結）で説明できる。根1（|Stab|=24）・根11（|Stab|=12）とも**非単位元すべてで非連結**。層ごとの分布 {1,1,1,1,1,1}（全層に固定セルがあるのに2成分）という 16胞体の鏡映と同型のケースも現れる。**定理は3種について述べられる。**

> ### ⚠️ 方法論上の落とし穴：セル隣接グラフの Aut は幾何的対称群より大きいことがある
>
> `o3x3x3o`（bitruncated 5-cell，10セル）は **\|Aut(セル隣接グラフ)\| = 3840 に対し多胞体の等長変換は240**。グラフを鵜呑みにすると上限 **48** が出るが，**正しくは 6**。
>
> 対策：安定化群を**セル重心への直交写像で実現される置換だけに制限**してから探索する（`isometryQ`：$A = \mathrm{LeastSquares}[C, C[[h]]]$ が直交かつ $C\!\cdot\!A = C[[h]]$）。この修正後，**全64種で極大 $H$ がすべて等長変換で実現されている**ことを確認済み。影響を受けたのは `o3x3x3o` 1種のみだが，論文では明記すべき。

#### 次の課題

- 24胞体・120胞体で「上限＝面安定化群」を**証明**する（8胞体は済み）
- 120胞体の最良木から STL を出して見た目を確認（未着手）
- 一様多胞体への適用：MIXED 2種・ALL INVALID 1種で対称な木が valid になるか
- 外接球半径が木に依らない理由（24/120胞体）

### 残る留保と未着手

- 3つの結果を分ける不変量の探索（§5 Open Problems に記載）
- MIXED 2種と x3o4o3x の自己交差構造の分析（`src/peel4d/analyze_120cell_overlap_bands.m` 相当）
- JCDCG³ 発表内容との整合（ユーザー方針：`And more results...` の形で部分的に新情報を入れる）
- 120/600胞体は BFS 木が 10⁴⁴・10⁶⁵ で全数列挙が原理的に不可能。**構造的証明が必要**。局所性（交差は木の距離が有界なセル間でしか起きない）が示せれば突破口になりうる

---

## 重要な幾何学的発見

### C2 候補の w 等値定理（2026-05-11）

任意の正多胞体の任意の C1 について，face-up（cell-centroid-up）後の C2 候補全員が同じ w 座標を持つ
（spread < 10⁻¹⁵）。理由：C1 の安定化部分群が C2 候補を推移的に置換し +w 軸を保存。

**含意**：RZ の "max w" 基準は k=2 で常にタイ → C2 選択は k=2 タイブレーカー（xy クロス積）のみで決まる。

### 600胞体は Impossible（構造的限界，2026-05-12 結論）

試したアプローチ：RZ/RS，fallback 有無，1段 BT，3D-face-centroid-up，vertex-up
→ **すべて 0/2400**

- cell-centroid-up：停止ステップ 5 固定（146, 150, 276, 279, 284）— 正二十面体対称のボトルネック
- 3D-face-centroid-up：停止ステップ 25 通りに多様化するが impossible 不変
- vertex-up：停止ステップ 14 通り，多くがより早く詰まる
- 停止ステップは cell-centroid-up 層境界に対応（276-284 ≈ 累積 273 = 赤道帯遷移）
- 根本原因：接続数 4（120胞体は 12）という貧しいグラフ構造 + Det 条件が全候補を封鎖

**結論**：グリーディ系では構造的に打破不可能。多段 BT は計算コスト的に非現実的。論文では「接続数 4 に由来する構造的限界」として説明済み。

### 120胞体 RZ nets の幾何学的多様性（C1=3）

12 通りの C2 すべて成功 → `netGeoKey4D` で **7 種類の geo-unique net**：
- 軌道：{1,15,17,97}（4），{4,18,51}（3），singleton 5 個（2, 10, 19, 100, 115）
- 普遍的内部コア：k=1..13 の r_3D が 7 net 全て一致，k=14 から分岐
- C2 重心は黄金比座標で**正二十面体**を形成（C1 が正十二面体セルの帰結）
- 6 対の対蹠ペア；対蹠は k=2 タイブレーカー符号反転で同軌道に属さない
- 螺旋パターン：型 A（CCW 上昇，5 net），型 B（CW 反転，C2=2），型 C（柱状，orbit {4,18,51}）

### 120胞体の自己交差構造（C1=3，全 12 C2 × 308 重なりペア）

- 同帯内 88 ペア（28.6%），異帯間 220 ペア（71.4%）
- 異帯重なりはほぼすべて隣接帯間
- 帯境界での代替セルは 1 個しか存在せず，グリーディ局所変更では回避不可
- 解消にはグローバル最適化が必要（グリーディ系では困難）

---

## 論文の現状（`paper_draft.tex`，2026-05-16 時点）

### 投稿準備状況

- **arXiv 投稿は即可**
- **CGTA（第1候補）は submittable，ただし major revision を覚悟すべき水準**
- **DCG レベルを狙うなら追加の形式化（命題化・600胞体の structural argument）が必要**

詳細は `commentsOpus260516.md` および `HISTORY.md` の「論文修正履歴」を参照。

### 査読で指摘されうる主要点（Opus 4.7 レビュー）

1. ~~等変性結果が「Remark」止まり~~ → **完了（2026-05-17）**：Remark 3.1 を `Proposition 3.1 [Equivariance and the 0/100% dichotomy]` + `Proof` に格上げし，実装詳細は残余 Remark 3.1（label `rem:symmetry3d` 維持）に整理。4 箇所のクロス参照を更新。`prop:equivariance` を新ラベルとして導入
2. ~~600-cell の "icosahedral bottleneck" が経験則止まり~~ → **完了（2026-05-17）**：Discussion「The 600-Cell」段落に Worked Example（`src/peel4d/analyze_600cell_stuck_example.m` で抽出した (C1=1, C2=2, k=146, last=cell 20) を表 `tab:600stuck-example` として）を追加。スタックの真の原因が「Det フィルタによる候補排除」ではなく「貪欲 max-w が外殻 4-隣接を先に消費し，唯一の下方出口も並行ブランチで既訪となるデッドエンド」であることを 4-正則性と関連付けて記述。120-cell（12-regular）との比較で構造的差を明示。Proposition 3.1 を介して 5 通り停止ステップへの propagation を説明
3. ~~コード/データ可用性ステートメントの欠如~~ → **完了（2026-05-17）**：`paper_draft.tex` の Acknowledgements 直後に `\section*{Code and data availability}` を追加；GitHub: <https://github.com/takashi-randomwalker/apple-peel-4d>
4. ~~参考文献が 10 件と寡少~~ → **完了（2026-05-17）**：4 件追加して 14 件に（`Pak2010` book draft, `Bern2003` Comput. Geom. 24 51-62, `Schlickenrieder1997` TU Berlin Diplomarbeit, `Coxeter1973` Dover 3rd ed.）。Towle は権威ある同名参照が存在しないためスキップ（Devadoss2022 が代替として既出）
5. ~~「なぜこの 2 規則か」の動機付けが弱い~~ → **完了（2026-05-17）**：Section 2.2 規則定義の直前に 1 段落追加。nearest neighbour / min angular deviation / smallest dihedral angle といった local rule が等変性を破る点を Proposition 3.1 に紐づけて説明し，global +z 軸と c_1 参照を使う rule の中で max azimuthal turn (RS) と max axial conservation (RZ) が 2 つの自然な端点であることを動機として記述

### Minor 改善項目（2026-05-17 セッション後半）

- ~~Abstract に Perfect/Possible/Impossible 分類が貢献の一部であることを明記~~ → **完了**：「A principal contribution is a three-way classification...」の文と「equivariance argument showing that face-transitive solids are confined to the 0/100% dichotomy」の補足を追加
- ~~Conclusion の Future Work をより具体的に~~ → **完了**：3 項目（two-stage, valid-net 特徴付け, random convex polyhedra）を 4 項目に拡張・具体化。新規追加：(1) hybrid rule の λ-parametrisation，(2) Darboux torsion / per-cell handedness による valid net 予測，(3) S²/S³ ランダム凸包での scaling 解析，(4) 600-cell の rigorous structural lemma（DCG レベル向け）
- ~~Section 5 重複整理~~ → **完了**：旧 5.2 (5-Cell) と 5.3 (16-Cell) を統合し「5-Cell and 16-Cell: rule agreement and partial coverage」（約 12 行）に圧縮。Section 4 / Table 4dglobal への参照に置換。5.1 (Cross-Dim Summary) と 5.4 (120-Cell) は独自内容のため保持。32 ページ維持
- ~~Acknowledgements ORCID / COI 宣言~~ → **完了**：Acknowledgements 末尾に `\paragraph{ORCID.}` で T. Yoshino: `0000-0003-1756-0162` を追加。標準 COI 宣言（"no known competing financial interests..."）も追加。S. Chaidee の ORCID は未取得（次回本人に確認予定）

### Sonnet 4.6 セカンドオピニオン対応（2026-05-17 セッション最終）

独立に走らせた Sonnet 4.6 レビューが新規 9 点を指摘 → 全て対応:

| 指摘 | 対応 |
|------|------|
| MW1: 3D Prop の RZ proof で `A(+z)=+z` 未トレース | 3 ステップ chain（+z → c_{F_1} → c_{σF_1} → +z）を proof に追加 |
| MW2: 600-cell 段落で Prop 3.1 が 5-step pattern を「説明する」と書いた overreach | "empirically" 表現に修正．orbit サイズ（742/825/404/323/106）と termination-step counts の対応を事実として記述．構造的導出は未達と明示 |
| MW3: 4D 等変性が一文 remark | `Proposition 4.1 [SO(4)-invariance]` + Proof + `Remark 4.1 [Partial equivariance and the role of the xy-plane]` を新設．label `prop:equivariance4d`，`rem:partialequivar4d` |
| MW4: 文献 entry 不完全 | Akitaya2024（著者 Samanta + Akitaya，FWCG24 URL），Devadoss2022（vol 111, article 101977, 2023, DOI），Kaino2019（11th Symmetry Congress Kanazawa の proc であることを追加） |
| MW5: k≥3 Det=0 fallback 動機なし | xy-area が「+w の orthogonal complement での azimuthal turn 指標」であること，4-点が共平面のとき volume form score が discriminate 不可になること，face-index tiebreak の非等変性を避けるための置換であることを 1 段落で説明 |
| MI3: Figure 5 caption の 16-cell uniqueness | 5/8/24-cell に加えて 16-cell も 20 successful pairs すべてが congruent であることを明示 |
| MI4: 「many square faces」claim が誤り | 7 つの 0% 解の真の共通項は **hexagonal faces 無し**（squares ではない）に修正．Icosidodecahedron 等が反例であることに対処 |
| MI6: 3D RS = min vs 4D RS = max の符号反転説明 | 「filter が admit する半空間の最も extreme value を選ぶ」共通原理として書き直し |
| MI8: §5.3 (120-cell) 冒頭が §4.2 と重複 | 冒頭 2 段落を 1 短段落に圧縮，Table 4dglobal への参照に置換．band structure / geo-diversity / self-intersection / face-centroid-up は §5.3 独自内容として保持 |

ページ数: 32 → 34（MW3 の Prop+Proof+Remark 追加で +1，MW5 の fallback motivation で +1）。

### Kaino2019 ページ番号

書誌調査エージェントが「あなたの draft = 142-145 vs arXiv:2604.16204 = 25-30」の conflict を検出したが，**著者本人が 142-145 を確認済み**（2026-05-19）．

### 2026-05-19 追加修正（Opus 4.7 第2回精査）

| 修正 | 内容 |
|------|------|
| Section 4.1 cell-centroid-up | 「face に整合しない・w 軸周りの残留自由度がある」を明記；3D face-up との比較（SO(2) 残留は不変）；(C1,C2) ペアが参照フレームを暗黙に固定する点を説明 |
| Section 6.1 Discussion | first-neighbor shell 段落を新設：全正多胞体で隣接セルの w が等値（対称性安定化部分群）→ max-w は識別不能 → geo-score が決定；4D Det が c1=(0,0,0,w1) の下で -w1·det_xyz に還元；120-cell の "universal inner core" を説明 |
| Table 7 本文言及 | Figure 4d-nets の後に Table summary への参照を自然な文脈で追加 |
| クロスリファレンス修正 | Section 5.1・Conclusion の `sec:4d` → `sec:discussion`（600-cell worked example の実際の位置） |
| US 綴り統一 | neighbour→neighbor（7箇所），realisation→realization（8箇所），analysed/summarised/characterising/minimising/maximising → 米英対応 |
| Section 4.3 | 「cells are large」→ 曲率蓄積による自己交差の正確な説明に修正 |
| Section 5.2 | "Across all 12 orderings" → "For a fixed starting cell" |
| "second band" 残留 | "second band onwards" → "once the algorithm leaves the first-neighbor shell (at k=14, ...)" |
| compare_3d4d_band1.m | 120-cell first-neighbor shell 検証スクリプト新規作成・実行 |

ページ数：34 → 35（first-neighbor shell 段落追加）．

### 推奨フロー

**現状（2026-05-19）**：Opus 4.7 第1回 5 点 + Minor 4 件 + Sonnet 4.6 追加 9 点 + Opus 4.7 第2回 9 点 = **全 27 点**を完了．Kaino ページ番号確認済み．DCG は狙わない方針．
→ **CGTA submittable**．arXiv 投稿はいつでも可能．S. Chaidee の ORCID が判明次第追記．

### 投稿候補ジャーナル

| 優先度 | ジャーナル |
|:------:|-----------|
| 第1候補 | Computational Geometry: Theory and Applications（Elsevier） |
| 第2候補 | Discrete & Computational Geometry（Springer） |
| 第3候補 | Journal of Computational Geometry |
| 代替 | Graphs and Combinatorics（Springer） |

### 論文構成（2026-05-16 時点）

```
1. Introduction
2. Algorithm（Geometric Setup / Two Selection Rules / Net Construction）
3. Results: 3D Polyhedra（Platonic / Archimedean / Mirror / TruncIcosa / SnubCube）
4. Extension to Four Dimensions（4D Algorithm / Results / 3D Realisation）
5. Computational Examples（Cross-Dim Summary / 5-cell / 16-cell / 120-cell）
6. Discussion（Face-Type Uniformity / Spiral vs Zonal / 600-cell / Companion）
7. Conclusion
```

### 注意事項

- bibitem キー：`\bibitem{Yoshino2026arXiv}`（arXiv 2604.16204）
- 3D RS フォールバック（Darboux Frame）は Algorithm 節で説明，4D アルゴリズムには使用しない
- xy-projection アプローチ（`src/peel4d/peeling4Df.m`）は **論文から削除済み**

---

## 可視化スクリプトの使い方

すべて Mathematica の `Get[...]` で読み込む。

### 3D 正多面体

```mathematica
Get["/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/visualize_platonic_nets_v2.m"]
```

`data/dataPlatonic.mx`（標準）と `data/dataPlatonic_strict.mx`（厳密）を読み込み，5種 × RS/RZ × withFallback/noFallback を描画。R2 は廃止済みのため非表示。

### 3D アルキメデス

```mathematica
Get["/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/visualize_archimedean_nets.m"]
```

`data/archimedean_faceup_results.mx` を読み込み，13 種 × RS/RZ × with fallback を描画。`showR2 = True` で R2 表示，`showR2 = False` がデフォルト。SnubCube は mirror（右手系）も inline 表示。

### 4D 正多胞体

```mathematica
Get["/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/unfold4D.m"]
```

`data/ans4DGlobal_v4.mx` + `src/peel4d/unfold3DExport.m` で 5/8/16/24/120 胞体の 4D→3D 展開図を `Graphics3D` 描画。等変性により C1=3 を代表として固定し全成功 C2 バリアントを表示。色付け：青（C1）→緑→赤（末尾）グラデーション。

### 設定パラメータ（スクリプト冒頭で変更可）

| パラメータ | デフォルト | 説明 |
|-----------|:----------:|------|
| `netSize` | 120 | 各展開図の画像サイズ (px) |
| `maxShowOK` / `maxShowFail` | `Infinity` | 表示上限 |
| `deduplicate` | `True` | 同じ order を 1 件に集約 |
| `deduplicateGeo` | `False` | 幾何学的重複除去（鏡像・回転等価を同一視）。ValueQ 保護済み |
| `showStrict` | `True` | 厳密条件結果も表示（Platonic 版） |
| `showR2` | `False` | R2 も表示（Archimedean 版） |

### `p3lUnfoldNet` の表示規約

- **多面体の外側から見たとき，螺旋が時計回り（CW）に見える**（右利きが左方向に剥く視点）
- 右半空間条件（det ≦ eps）で RS は 3D 空間で外側から CW の経路を直接選ぶため，**表示での y 反転は不要**
- 4D 展開図にはキラリティ正規化を適用しない（Graphics3D は視点自由のため不要）

---

## コードを読むときの注意

- `summary.tex`（または `summary_en.tex`）が最も包括的な情報源
- 3D版の右条件は `Det[{cc[[top]], cc[[last]], cc[[j]]}] <= eps`（c_1 大域固定参照，`eps = 10^-10`）。`src/peel3d/peeling3DLoxo.m` の `p3lPeelPair` 参照
- 3D版の選択基準も Det ベース（2026-05-12）：RS は min Det，RZ は max z → min Det タイブレーク。RS フォールバックのみ Darboux frame min φ
- **単一候補優先は廃止**（2026-05-06）：大域螺旋保証のため，単一候補でも右条件を適用
- 4D版（`src/peel4d/peeling4Df4.m`）は k=2 で xy クロス積，k≥3 は Det，Det=0 縮退時は xy クロス積にフォールバック
- 旧版 `src/peel4d/peeling4Df.m` / `src/peel4d/peeling4Df3.m` は左条件が異なる（xy 2成分 / xyz ローカル参照）— **使用非推奨**
- アルキメデス計算結果は 2026-05-12 に Det ベース選択で再計算済み
- Mathematica の `Round[x, N]` は N 小数点以下ではなく N の最近傍倍数に丸める（3 桁表示には `Round[x, 0.001]`）

---

## コンパニオン論文（arXiv:2604.16204）との実装差異

`src/peel3d/peeling3DLoxo.m`（本論文）は旧コード `peeling3Df` と次の 3 点で異なる。これらは等変性の回復を目的とした修正。

| 項目 | 旧 `peeling3Df` | 現行 `src/peel3d/peeling3DLoxo.m` |
|------|----------------|----------------------|
| フィルタ参照点 | c_k（局所，毎ステップ変化），厳密 `> 0` | c_1（大域固定），`<= ε`，ε = 10⁻¹⁰ |
| タイブレーク | リスト順（`Position` の最初の要素，非幾何） | min Det（幾何学的・等変） |
| RS 規則 | なし | あり（min Det 選択） |

旧コードは正十二面体で 53.3%（等変性理論は 0% か 100% のみ）という矛盾を生んでいた。
詳細な経緯は `HISTORY.md` の「数値誤差問題の解決」節を参照。
