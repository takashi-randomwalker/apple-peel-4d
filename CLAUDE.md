# CLAUDE.md — 4D Apple-Peel 展開プロジェクト

## プロジェクト概要

3次元多面体の **Apple-Peel 展開**（Yoshino, Chaidee）を **4次元正多胞体**に拡張するプロジェクト。

- 実装言語：**Mathematica**（`.m` スクリプト + `.nb` Notebook）
- ドキュメント：LaTeX / Markdown
- 対象：6種の正多胞体（5胞体・8胞体・16胞体・24胞体・120胞体・600胞体）

詳細な計算ログ・廃止された実装の比較・セッション別修正履歴は **`HISTORY.md`** を参照。
アルゴリズムの理論的詳細・等変性証明は **`summary.tex` / `summary_en.tex`** を参照。

---

## アルゴリズム（現行版）

### 3次元版（`peeling3DLoxo.m`）

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
`run_platonic_updated.m` の `p3lPeelPairNoFB2`：ステップ 4 でフォールバック使わず即終了（右条件は適用）。

### 4次元版（`peeling4Df4.m`，現行・推奨）

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

> **旧版 `peeling4Df.m` / `peeling4Df3.m`** は左条件が異なる（xy 2成分 / xyz ローカル参照）。現役は `peeling4Df4.m` のみ。論文からも旧版の議論は削除済み。

#### 1段バックトラック（`peeling4Df4_1back.m`）

詰まったとき1ステップ戻って別候補を試す拡張（justBT フラグで連続バックトラック禁止）。全6胞体で **グリーディ版と完全に同一結果**（改善なし，2026-05-11）。

---

## 主要ファイル（現役）

### 4D 本体・実行スクリプト

| ファイル | 役割 |
|----------|------|
| `peeling4Df4.m` | 4D本体（c1–c2 平面グローバル参照 + RZ/RS） |
| `peeling4Df4_1back.m` | 1段バックトラック版（グリーディと同値） |
| `run4DGlobal_update.m` | 5/8/16/24胞体の v4 計算（k≥3 Det=0 fix 適用） |
| `run4DGlobal_120cell.m` | 120胞体の v4 計算 |
| `run4DGlobal_nofallback.m` | RZ/RS × with/without fallback 比較 |
| `run4DGlobal_RSRZ.m` | RS vs RZ 比較 |
| `run4DBacktrack.m` | 1段バックトラック計算 |
| `run4DExportSTLv4.m` | v4 から valid net を STLs_v4/ に出力 |

### 3D 本体・実行スクリプト

| ファイル | 役割 |
|----------|------|
| `peeling3DLoxo.m` | 3D本体（RS/RZ + R2） |
| `run_platonic_updated.m` | 正多面体 5種計算（fallback 有無） |
| `run_platonic_strict.m` | 厳密条件（eps<0）計算 |
| `run_archimedean_faceup.m` | アルキメデス 13種（fallback あり） |
| `run_archimedean_nofallback.m` | アルキメデス 13種（fallback なし） |

### 可視化・展開・STL 出力

| ファイル | 役割 |
|----------|------|
| `unfold3DExport.m` | SVD射影・Procrustes整列・SAT判定・STL出力 |
| `unfold4D.m` | 4D展開図の3D可視化（`ans4DGlobal_v4.mx` + `unfold3DExport.m`） |
| `visualize_platonic_nets_v2.m` | 正多面体展開図描画 |
| `visualize_archimedean_nets.m` | アルキメデス多面体展開図描画 |
| `gen_120cell_faceup_nets.m` | 120胞体 cell-up vs face-up 比較（v4 fix 適用） |

### データファイル（現役）

| ファイル | 内容 |
|----------|------|
| `ans4DGlobal_v4.mx` | **現行最新**：4D 全6胞体の計算結果（v4） |
| `ans4DGlobal_nofb.mx` | フォールバックなし比較（サマリのみ） |
| `ans4DGlobal1back.mx` | 1段バックトラック結果 |
| `dataPlatonic.mx` | 正多面体 × R1/R2/R3 × fallback 有無 |
| `dataPlatonic_strict.mx` | 正多面体 × RS/RZ × fallback 有無（eps<0） |
| `archimedean_faceup_results.mx` | アルキメデス（fallback あり，Det ベース） |
| `archimedean_nofallback_results.mx` | アルキメデス（fallback なし，Det ベース） |
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

> 旧版（`peeling4Df.m`, `peeling4Df3.m`, `peeling1back.m`, `ans4DGlobal.mx`, `ans4DGlobal_v2.mx`, `ans4DGlobal_v3.mx`, `STLs/`, `STLs_v3/` など）は `HISTORY.md` 参照。

---

## 計算結果（最新確定）

### 4D：`peeling4Df4.m`（v4，k≥3 Det=0 fix，2026-05-12）

結果：`ans4DGlobal_v4.mx`

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

### フォールバックなし比較（`peeling4Df4.m` noFB）

| 多胞体 | RZ(w) | RZ(n) | RS(w) | RS(n) |
|--------|------:|------:|------:|------:|
| 5/8/24胞体 | Perfect | Perfect | Perfect | Perfect |
| 16胞体 | Possible/20 | Possible/20 | Impossible | Impossible |
| 120胞体| **Perfect/1440** | Possible/978 | Impossible | Impossible |
| 600胞体| Impossible | Impossible | Impossible | Impossible |

- フォールバックが有効な唯一の組み合わせは **120胞体 × RZ**（462ペア=32%がフォールバック経由）
- 4次元では RZ が RS を大幅に上回る（RS は 16胞体・120胞体で Impossible）

### 3D 正多面体（`dataPlatonic.mx`，2026-05-06）

全結果が 0% か 100%（面推移性による等変性理論と一致）。
- with fallback：全 5 種・全規則で 100%
- no fallback：Tetrahedron のみ全規則で 100%；他は R2 のみ 0%（他規則は 100%）
- 厳密条件（eps<0）：Tetrahedron のみ no-fallback で 100%，他は 0%

### 3D アルキメデス（`archimedean_faceup_results.mx`，2026-05-12 Det ベース）

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
| `run_random_polyhedra.m` | N=30 の Delaunay vs Voronoi 詳細比較（1 インスタンス + 20 試行） |
| `run_random_scaling.m` | N = {20,30,50,100,200} × 4 条件のスケーリング実験 |
| `run_repulsive_comparison.m` | Thomson Coulomb 緩和による均一分布の効果実験（N=30） |

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

### 計算結果（`face_rotation_net_all4D_v2.m`，2026-06-09）

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
| `face_rotation_net_all4D_v2.m` | 全6胞体の BFS ネット検証（体積基準，**現行推奨**） |
| `face_rotation_net_small4D.m` | 5/8/16/24-cell 検証（v1，RegionMember 版，歴史的参考） |
| `face_rotation_net_600cell.m` | 600-cell 専用検証スクリプト |
| `check_16cell_detail.m` | バグ調査：RegionMember vs 体積基準の比較 |
| `face_rotation_net_120cell.m` | 120-cell BFS ネット可視化（net PNG 出力） |
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
| Open Problems 段落削除 | 「apple-peel 自己交差は distant layer」記述は実測（隣接帯間が主体，`analyze_120cell_overlap_bands.m` 由来）と矛盾するため段落ごと削除 |
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
| `uniform4D_wythoff.m` | **生成器**：Coxeter 図から一様4-多胞体の頂点・辺・面・セルを生成（A₄/B₄/F₄，33種） |
| `run_wythoff_bfs.m` | 生成データを BFS チェッカに流す（Phase A：正多胞体を全根／Phase B：全種を root 1） |
| `run_wythoff_bfs_allroots.m` | **全30種 × 全根**（1893根）の確定計算。結果を `wythoff_bfs_allroots.mx` に保存 |
| `check_runcinated24.m` | 唯一の例外 x3o4o3x の全240根詳細（重なりペア数の分布） |

`face_rotation_net_all4D_v2.m` は本体の main ループを走らせずに定義だけ読み込む（テキストを `Print["Face-rotation BFS net check v2` の直前で切って `ToExpression`）。

### 生成アルゴリズム（`uniform4D_wythoff.m`）

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

`face_rotation_net_all4D_v2.m` が使うのは [[1]], [[4]], [[6]] のみ。

- [[3]] を「セルの頂点リスト」と誤読しやすい：f5.m では頂点数5・次数4がセル数5・4頂点と**偶然一致**する
- **`f8.m` の [[2]] は32本の辺を両方向で64件**列挙（f5/f16/f24 は1回ずつ）。[[2]] は誰も読まないので実害なし。生成側は無向1回の規約

### H₄ 系（2026-08-12 追加）

ridge 検出を **O(C²) → O(|W|·多角形サイズ)** に書き換えて H₄（位数14400）に到達。

**書き換えの中身**（`uniform4D_wythoff.m`）
- **ridge**：Wythoff 的多胞体の d-面は階数 d の放物型部分群の軌道。rank-2 放物型 W_{i,j} で ridge の原型を作り，その**頂点集合**を W で軌道に乗せる。セル対の総当たりが不要に
- **ridge が属する2セル**：頂点→セルの逆引きを作り，ridge の全頂点が属するセルの共通部分（ちょうど2個になるはず）
- **セル所属判定**：超平面ごとに1回の行列ベクトル積にベクトル化
- 検証項目が2つ増加：**各 ridge がちょうど2セルに属する**，**全頂点の次数が一様（頂点推移性）**
- A₄/B₄/F₄ の33種は書き換え前と**完全に同一の V/E/F/C**（回帰テスト通過）かつ高速化

**生成結果：48/48 が全自己検証を通過**（A₄ 9 + B₄ 15 + F₄ 9 + H₄ 15），重複3組を除いて**正味45種**。総生成時間509秒，`wythoff_gen_all.mx` に保存。

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

**MIXED の裏取り（`verify_mixed_h4.m`，2026-08-13）**：数値誤差の疑いを排除済み。

- 失敗根5個・成功根3個 × 2多胞体 = 16根で **SAT と RegionMeasure が判定・重なり数とも完全一致**
- eps を 10⁻⁴ 〜 10⁻¹⁰ まで振っても重なり数は不変（閾値上の縁の判定ではない）
- x5o3x3x の重なり数分布は滑らか：`{{0,47},{1,282},{2,361},{3,546},{4,608},{5,439},{6,225},{7,76},{8,34},{9,12},{10,9},{11,1}}`
- x5x3o3x は `{{0,2570},{1,70}}` で，失敗根はすべて重なり1個ちょうど

### snub 24-cell と grand antiprism（2026-08-12，`uniform4D_special.m`）

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

**BFS 全根：両方とも ALL VALID，重なり0**（144/144 が31.6秒，320/320 が344.9秒）。結果は `special_bfs_allroots.mx`。

### 角柱17種（2026-08-13，`uniform4D_prisms.m`）

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

`export_uniform_data.m` で64種のデータを出力。**64ファイル 13.3 MB**。`index.csv` 付き。

**5要素形式**（`_4DData/f*.m` の6要素から2つ削り1つ足した）：

```
[[1]] verts        頂点座標（機械精度）
[[2]] edges        頂点添字ペア
[[3]] faces        2-面（頂点添字，巡回順）
[[4]] cellsByFace  セル（面添字）
[[5]] cellAdj      各セルの隣接セル  ← 新規
```

**削った2つ**：`vertAdj`（第3要素）と `faceAdj`（第5要素）。理由は容量ではない。
- **`face_rotation_net_all4D_v2.m` は `raw[[1]], raw[[4]], raw[[6]]` の3つしか読まない**（138行目）。この2つはリポジトリ内のどこからも参照されていない
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

`export_uniform_nets.m` で **67ケース**（64種の root 1 ＋ 例外3種の追加根）を出力。合計55MB。

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

**Figure 2（`fig:uniformfail`）新設**：`make_uniform_figure.m` で生成。失敗の2つの型を1つずつ。
- (a) runcinated 24-cell の**最良の根**（root 47，交差3組）— 「どの根でもダメ」は最悪例より最良例が効く
- (b) runcitruncated 120-cell の失敗根（root 748）の**交差ペア近傍36セルの拡大図**— 2640セル全体は団子になるため
- 交差ペアの2セルは**赤と青に塗り分け**（同色だと相貫が1個の多面体に見える）

Abstract に1文追加（64種24,487根，61/1/2）。Introduction の構成説明に §4 を追加。Code and Data Availability に角柱・SAT・`uniform_nets/` を追記。**9ページ**でビルド通過，警告なし。

### 残る留保と未着手

- **タイブレーク規則依存**：BFS 木はタイブレーク（現行はセル添字順）に依存。セル推移的なら無関係だが一様多胞体は複数のセル軌道を持つため，MIXED の根の個数（47/2640 等）は**この規則の下での値**。論文 §4 Caveats に明記済み
- 3つの結果を分ける不変量の探索（§5 Open Problems に記載）
- MIXED 2種と x3o4o3x の自己交差構造の分析（`analyze_120cell_overlap_bands.m` 相当）
- JCDCG³ 発表内容との整合（ユーザー方針：`And more results...` の形で部分的に新情報を入れる。アウトライン確定後に議論）

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
2. ~~600-cell の "icosahedral bottleneck" が経験則止まり~~ → **完了（2026-05-17）**：Discussion「The 600-Cell」段落に Worked Example（`analyze_600cell_stuck_example.m` で抽出した (C1=1, C2=2, k=146, last=cell 20) を表 `tab:600stuck-example` として）を追加。スタックの真の原因が「Det フィルタによる候補排除」ではなく「貪欲 max-w が外殻 4-隣接を先に消費し，唯一の下方出口も並行ブランチで既訪となるデッドエンド」であることを 4-正則性と関連付けて記述。120-cell（12-regular）との比較で構造的差を明示。Proposition 3.1 を介して 5 通り停止ステップへの propagation を説明
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
- xy-projection アプローチ（`peeling4Df.m`）は **論文から削除済み**

---

## 可視化スクリプトの使い方

すべて Mathematica の `Get[...]` で読み込む。

### 3D 正多面体

```mathematica
Get["/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/visualize_platonic_nets_v2.m"]
```

`dataPlatonic.mx`（標準）と `dataPlatonic_strict.mx`（厳密）を読み込み，5種 × RS/RZ × withFallback/noFallback を描画。R2 は廃止済みのため非表示。

### 3D アルキメデス

```mathematica
Get["/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/visualize_archimedean_nets.m"]
```

`archimedean_faceup_results.mx` を読み込み，13 種 × RS/RZ × with fallback を描画。`showR2 = True` で R2 表示，`showR2 = False` がデフォルト。SnubCube は mirror（右手系）も inline 表示。

### 4D 正多胞体

```mathematica
Get["/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/unfold4D.m"]
```

`ans4DGlobal_v4.mx` + `unfold3DExport.m` で 5/8/16/24/120 胞体の 4D→3D 展開図を `Graphics3D` 描画。等変性により C1=3 を代表として固定し全成功 C2 バリアントを表示。色付け：青（C1）→緑→赤（末尾）グラデーション。

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
- 3D版の右条件は `Det[{cc[[top]], cc[[last]], cc[[j]]}] <= eps`（c_1 大域固定参照，`eps = 10^-10`）。`peeling3DLoxo.m` の `p3lPeelPair` 参照
- 3D版の選択基準も Det ベース（2026-05-12）：RS は min Det，RZ は max z → min Det タイブレーク。RS フォールバックのみ Darboux frame min φ
- **単一候補優先は廃止**（2026-05-06）：大域螺旋保証のため，単一候補でも右条件を適用
- 4D版（`peeling4Df4.m`）は k=2 で xy クロス積，k≥3 は Det，Det=0 縮退時は xy クロス積にフォールバック
- 旧版 `peeling4Df.m` / `peeling4Df3.m` は左条件が異なる（xy 2成分 / xyz ローカル参照）— **使用非推奨**
- アルキメデス計算結果は 2026-05-12 に Det ベース選択で再計算済み
- Mathematica の `Round[x, N]` は N 小数点以下ではなく N の最近傍倍数に丸める（3 桁表示には `Round[x, 0.001]`）

---

## コンパニオン論文（arXiv:2604.16204）との実装差異

`peeling3DLoxo.m`（本論文）は旧コード `peeling3Df` と次の 3 点で異なる。これらは等変性の回復を目的とした修正。

| 項目 | 旧 `peeling3Df` | 現行 `peeling3DLoxo.m` |
|------|----------------|----------------------|
| フィルタ参照点 | c_k（局所，毎ステップ変化），厳密 `> 0` | c_1（大域固定），`<= ε`，ε = 10⁻¹⁰ |
| タイブレーク | リスト順（`Position` の最初の要素，非幾何） | min Det（幾何学的・等変） |
| RS 規則 | なし | あり（min Det 選択） |

旧コードは正十二面体で 53.3%（等変性理論は 0% か 100% のみ）という矛盾を生んでいた。
詳細な経緯は `HISTORY.md` の「数値誤差問題の解決」節を参照。
