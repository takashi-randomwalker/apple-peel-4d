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
