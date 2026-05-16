# CLAUDE.md — 4D Apple-Peel 展開プロジェクト

## プロジェクト概要

3次元多面体の **Apple-Peel 展開**（Yoshino, Chaidee）を
**4次元正多胞体**に拡張するプロジェクト。

- 実装言語：**Mathematica**（`.m` スクリプト + `.nb` Notebook）
- ドキュメント：LaTeX / Markdown
- 対象：6種の正多胞体（5胞体・8胞体・16胞体・24胞体・120胞体・600胞体）

---

## アルゴリズムの要点

### 3次元版（`peeling3DLoxo.m`）

#### 前処理
- **Face-up 回転**：`p3lAlignTopToZ[vers, faces, top]` で開始面 F1 の重心を +z 軸に揃える
- **全ペア評価**：`p3lRunAll[vers, faces, rule]` が全 (F1, F2) ペアを評価

#### 1ステップの選択
1. `last` の未訪問隣接面集合を `pending` とする
2. `pending` が空なら終了（失敗）
3. **大域的右半空間条件**で候補を絞る（単一候補でも適用する）：
   - `Det[{c_1, c_k, c_j}] <= eps` を満たす面を右候補 `rightCands`（変数名は `leftCands`）とする
   - `c_1 = cc[[top]]`（開始面の重心，固定参照）；`eps = 10^-10`（数値誤差吸収）
   - 右半空間 = 外側から見て CW 方向の候補（右利きが左方向に剥く向き）
4. `rightCands` が空でなければ規則 r で次の面を選択；空なら **フォールバック**（pending 全体から規則 r の変形で選択）

**大域参照の理由**：c_{k-1}（局所参照）は各ステップで変化するため，螺旋方向の一貫性を保証できない。c_1（固定参照）を用いることで，全ステップで z 軸周りの大域的な螺旋方向を保証する。また，単一候補でも右条件を適用する（大域的螺旋保証のため）。

#### 選択規則の呼び名

> **⚠️ 論文方針（2026-05-10）**：`paper_draft.tex` では **R2 規準は掲載しない**。
> R2 に関する記述・列・例・サブセクションは論文に追加しないこと。
> R2 の計算結果や考察は記録として `summary.tex`・`summary_en.tex`・`CLAUDE.md` にのみ保持する。

論文で扱う2規則のみを議論するときは以下の呼び名を使う：

| 記号 | 正式名 | 略称 |
|------|--------|------|
| R1 | **Spiral rule** | **RS** |
| R3 | **Zonal rule** | **RZ** |

Spiral rule：最大右旋回（min Det）で外側から見て CW の密な螺旋軌跡を描く。
Zonal rule：最高緯度優先（max z）で緯度帯を一周してから降りる帯状走査。

#### 選択規則（2026-05-12 改訂：Det ベース統一）
| 規則 | 右候補あり | フォールバック（右候補なし） |
|------|-----------|---------------------------|
| R1 = RS（Spiral, min Det） | `rightCands` から min `Det[{c_1,c_k,c_j}]`（最も負 = 最大右旋回） | `pending` から Darboux frame の最小 φ |
| R2（loxodrome） | `rightCands` から max Det（最も 0 に近い = 最小右旋回） | `pending` から最大 \|Det\| |
| R3 = RZ（Zonal, max z + min Det 二次） | `rightCands` から最大 z；z 同値（差 ≦ 10⁻¹⁰）なら min Det で二次選択 | `pending` から最小 z；z 同値なら min Det で二次選択 |

**Det ベース選択の理由（2026-05-12）**：
`Det[{c_1,c_k,c_j}] = Cross[c_1,c_k]·c_j`（外積の内積）は右条件フィルタと同じ量であり，
フィルタと選択基準が完全に一貫する。局所 Darboux frame（c_{k-1}→c_k）は各ステップで変化するため
角度 φ が大域螺旋方向と必ずしも一致しないが，Det は c_1 固定参照なので全ステップで一貫する。
RS フォールバック（右候補なし時）のみ Darboux frame min φ を維持（フォールバック自体が非幾何的な救済）。

**RZ 二次基準の理由**：正十二面体のように高対称な多面体では，Face-up 後に複数の候補面が理論上同一 z を持つ。
浮動小数点誤差で z 値にわずかな差（~10⁻¹⁶）が生じリスト順でタイブレークされると等変でなくなる。
min Det を二次基準にすることで等変性を回復する（旧版 min φ と実質等価）。

- Darboux フレーム：`p3lDarbouxFrame[cPrev, cCurr]` → `{nHat, fHat, lHat}`（RS フォールバックのみ使用）
- 角度計算：`p3lAngleFromForward[fHat, lHat, cj]` で φ を計算（RS フォールバックのみ）

#### フォールバックなし版（`run_platonic_updated.m` の `p3lPeelPairNoFB2`）
- ステップ 4 でフォールバックを使わず即終了する変種
- 大域的右半空間条件（c_1 固定参照，eps）は適用する

### 4次元版（`peeling4Df.m` / `peeling4Df3.m` / `peeling4Df4.m`）

- **剥き軸**：w 軸（Face-up 配向）
- **左条件（現行 `peeling4Df.m`）**：`c_{k,y} * c'_x - c_{k,x} * c'_y >= 0`（xy 2成分外積）
- **左条件（`peeling4Df3.m`）**：`Det[{c_{k-1}, c_k, c_j}]_{xyz} >= 0`（xyz 3成分，ローカル k-1 参照，k≥2）
  - k=1 は直前セル未定義のため xy 2成分にフォールバック
- **左条件（`peeling4Df4.m`，推奨）**：`Det[{c_1, c_2, c_k, c_j}] >= -eps`（4次元体積形式，**c1–c2 平面グローバル参照**）
  - 3D 版 `Det[{c_1, c_k, c_j}] >= -eps` の直接類推
  - Face-up 後 c_1 = (0,0,0,w_1) で展開すると `-w_1 * Det_xyz[{c_2, c_k, c_j}] >= -eps` と等価
  - k=2 のとき Det[{c1,c2,c2,c_j}] = 0 → 全候補通過（特別扱い不要）
  - **等変性**：A ∈ SO(4) は det(A)=1 なので `Det[{Ac1,Ac2,Ac_k,Ac_j}] = Det[{c1,c2,c_k,c_j}]` が成立
- **選択規則**：規則3（左候補から max w；w 同点（差 ≦ 10^-10）なら max Det タイブレーク；フォールバック min w，同点なら max |Det|）
  - 3D RZ の max-φ 二次基準に対応する4Dタイブレーク
  - k=2 では c_k = c_2 のため Det[{c1,c2,c2,cj}] ≡ 0 → タイブレークに xy クロス積 `c2_x·cj_y − c2_y·cj_x` を使用
  - **k≥3 Det=0 フォールバック（2026-05-12）**：タイブレーク候補の Det が全て ≈0 のとき xy クロス積にフォールバック（k=2 と同じ式）。8胞体・120胞体でリスト順依存を解消
- **単一候補優先**：4D版には未実装（変更していない）

### 1段バックトラック（`peeling1back.m`）

- 詰まったとき1ステップ戻って別候補を試す拡張版

---

## 主要ファイル

| ファイル | 役割 | 状態 |
|----------|------|------|
| `peeling4Df.m` | 4D本体（xy 2成分左条件 + 規則3） | 完成・使用中 |
| `peeling4Df3.m` | 改良版（xyz ローカル k-1 参照 + 規則3） | 作成済み・**未テスト** |
| `peeling4Df4.m` | 4D版（c1–c2 平面グローバル参照 + 規則3）3D の素直な拡張 | 完成・使用中 |
| `run4DPeeling.m` | 全6胞体実行スクリプト（`peeling4Df.m` 用） | 完成 |
| `run4DPeelingR3.m` | `peeling4Df3.m` 用実行スクリプト | 作成済み・未実行 |
| `run4DPeelingGlobal.m` | `peeling4Df4.m` 用実行スクリプト（タイブレーク前版） | 完成・実行済み |
| `run4DFixed_top3.m` | C1=3 固定で全6胞体の隣接 C2 を計算，旧版との差分レポート | 完成・実行済み（2026-05-10） |
| `run4DGlobal_tiebreak.m` | タイブレーク版 peeling4Df4 を全 C1 で計算，`ans4DGlobal_v2.mx` 保存 | 完成・実行済み（2026-05-10） |
| `run4DGlobal_xy2.m` | xy クロス積 k=2 タイブレーク版（v3）で全6胞体・全C1計算 | 完成・実行済み（2026-05-10） |
| `run4DGlobal_update.m` | k≥3 Det=0 修正版で 5/8/16/24胞体を再計算，120胞体は旧データ流用 | 完成・実行済み（2026-05-12） |
| `run4DGlobal_120cell.m` | 120胞体を現行 peeling4Df4 で計算し ans4DGlobal_v4.mx に追記 | 完成・実行済み（2026-05-12） |
| `run4DGlobal_nofallback.m` | RZ/RS × with/without fallback で全6胞体を計算・比較 | 完成・実行済み（2026-05-12） |
| `run4DExportSTLv4.m` | v4 結果から valid net を STLs_v4/ に出力（重なり判定付き） | 完成・実行済み（2026-05-12） |
| `ans4DGlobal_nofb.mx` | フォールバックなし比較の計算結果（サマリのみ，フルデータなし） | 保存済み（2026-05-12） |
| `STLs_v4/` | 280枚の STL ファイル（valid net のみ，重なりなし保証） | 保存済み（2026-05-12） |
| `run4DGlobal_RSRZ.m` | RZ と RS の両規則で全6胞体・全C1計算・比較 | 完成・実行済み（2026-05-10） |
| `analyze_600cell_stuck.m` | 600胞体の全(C1,C2)ペアの終了ステップ頻度分布を解析 | 完成・実行済み（2026-05-10） |
| `analyze_600cell_wdist.m` | 600胞体の cell-centroid-up 後 w 分布と Kaino 層構造を比較 | 完成・実行済み（2026-05-12） |
| `test_600cell_vertexup.m` | 600胞体に vertex-up を適用（C1=1，16ケース）；全 STUCK を確認 | 完成・実行済み（2026-05-12） |
| `ans4DGlobal.mx` | `peeling4Df4.m` タイブレーク前の計算結果（Mathematica バイナリ） | 保存済み |
| `ans4DGlobal_v2.mx` | `peeling4Df4.m` タイブレーク後の計算結果（0バイト，DumpSave 失敗） | 注意：0バイト |
| `peeling1back.m` | 旧版1段バックトラック（xy ローカル参照，使用非推奨） | 参照用のみ |
| `peeling4Df4_1back.m` | 1段バックトラック版（Det グローバル参照 + RZ/RS，`peeling4Df4.m` 準拠） | 完成・実行済み（2026-05-11） |
| `run4DBacktrack.m` | 1段バックトラックで全6胞体・全C1計算，グリーディ版と比較 | 完成・実行済み（2026-05-11） |
| `ans4DGlobal1back.mx` | 1段バックトラック版の計算結果（Mathematica バイナリ） | 保存済み |
| `peeling3DLoxo.m` | 3D版 R2/R3 実装 | 完成 |
| `unfold3DExport.m` | SVD射影・Procrustes整列・SAT判定・STL出力（`unfoldTo3DFast`・lazy `satOverlap` 最適化済み） | 完成（2026-05-11性能改善） |
| `run4DExportSTLv3.m` | RZ全5胞体でペーリング再計算→SAT重なり判定→valid net STL出力 | 完成・実行済み（2026-05-11） |
| `unfold4D.m` | 4D展開図の3D可視化スクリプト（`ans4DGlobal_v4.mx` + `unfold3DExport.m` を使用）；C1=3 の全12 C2 を表示（v4） | 完成（2026-05-08，色・エッジ更新 2026-05-13） |
| `gen_120cell_faceup_nets.m` | 120胞体 C1=3 全12 C2 の cell-centroid-up vs 3D-face-centroid-up 比較（v4 fix 適用）；PDF + ノートブック inline 表示 + STL 出力 | 完成（2026-05-13） |
| `STLs_120cell_faceup/` | face-up 成功7ケース × cell-up/face-up の STL ファイル（14枚，重なりあり） | 保存済み（2026-05-13） |
| `ans4D.mx` | 現行版の計算結果（Mathematica バイナリ） | 保存済み |
| `dataPlatonic.mx` | 正多面体 5種 × R1/R2/R3 × fallback あり/なし の全ペア計算結果 | `run_platonic_updated.m` 実行後に生成 |
| `dataPlatonic_strict.mx` | 正多面体 5種 × RS/RZ × fallback あり/なし（eps < 0，thresh = +10^-10） | `run_platonic_strict.m` 実行後に生成 |
| `run_platonic_strict.m` | eps < 0（厳密条件）での RS・RZ 計算スクリプト | 完成・実行済み（2026-05-07） |
| `run_platonic_eps_comparison.m` | thresh の3ケース（eps>0/=0/<0）を比較するスクリプト | 完成・実行済み |
| `dataSummaryPlatonic.md` | 正多面体計算データの現状説明（保存形式・データ構造・可視化） | 作成済み |
| `visualize_platonic_nets.m` | 正多面体展開図の可視化スクリプト（`dataPlatonic.mx` を読み込んで表示） | 完成・使用中 |
| `run_archimedean_nofallback.m` | アルキメデス13種をフォールバックなし（RS/RZ）で計算するスクリプト | 完成・実行済み（2026-05-10） |
| `archimedean_nofallback_results.mx` | フォールバックなし計算結果（Mathematica バイナリ） | 保存済み |
| `visualize_archimedean_nets.m` | アルキメデス多面体展開図の可視化スクリプト（`archimedean_faceup_results.mx` を読み込んで表示） | 完成 |
| `STLs/` | 570枚の STL ファイル（旧版，重なり未チェック） | 保存済み |
| `STLs_v3/` | 285枚の STL ファイル（valid net のみ，重なりなし保証） | 保存済み（2026-05-11） |
| `_4DData/f5.m` 〜 `f600.m` | 正多胞体データ（頂点・面・セル） | 完成 |
| `summary.tex` / `.pdf` | アルゴリズム詳細ドキュメント（最も包括的） | 最新版 |
| `summary_en.tex` / `.pdf` | 英語版 summary（変数定義付き） | 作成済み |
| `archimedean_results.tex` / `.pdf` | Archimedean 13種での R1/R2/R3 比較結果 | 作成済み |
| `literature_review_spiral_unfolding.md` / `.tex` / `.pdf` | 螺旋展開の既往研究レビュー | 作成済み |
| `rule_comparison_illustration.pdf` | R1/R2/R3 の違いを示す正十二面体図 | 作成済み |
| `TODO_4D_resume.md` | 再開手順・詳細タスクリスト | 参照のこと |

---

## 計算結果

### `peeling4Df.m`（xy 2成分外積，旧版）

| 多胞体 | セル数 | 隣接数/セル | (C1,C2)総数 | ユニーク成功数 | 分類 |
|--------|:------:|:-----------:|:-----------:|:-------------:|------|
| 5胞体  | 5      | 4           | 20          | 20            | **Perfect** |
| 8胞体  | 8      | 6           | 48          | 48            | **Perfect** |
| 16胞体 | 16     | 4           | 64          | 26            | Possible |
| 24胞体 | 24     | 8           | 192         | 119           | Possible |
| 120胞体 | 120   | 12          | 1,440       | 357           | Possible |
| 600胞体 | 600   | 4           | 2,400       | —             | **未完了** |

### `peeling4Df4.m`（c1–c2 平面グローバル参照 + 幾何学的タイブレーク，2026-05-10）

計算スクリプト：`run4DGlobal_xy2.m`，結果：`ans4DGlobal_v3.mx`（参照用；現在は v4 が最新）

| 多胞体 | (C1,C2)総数 | True | False | Unique | 分類 |
|--------|:-----------:|-----:|------:|-------:|------|
| 5胞体  | 20          | 20   | 0     | 20     | **Perfect**   |
| 8胞体  | 48          | 48   | 0     | 48     | **Perfect**   |
| 16胞体 | 64          | 25   | 39    | 25     | Possible      |
| 24胞体 | 192         | 192  | 0     | 192    | **Perfect**   |
| 120胞体 | 1,440      | 1,250| 190   | 1,250  | Possible      |
| 600胞体 | 2,400      | 0    | 2,400 | 0      | **Impossible** |

### `peeling4Df4.m`（k≥3 Det=0 フォールバック修正版，2026-05-12）

計算スクリプト：`run4DGlobal_update.m`（5/8/16/24胞体）+ `run4DGlobal_120cell.m` + `run4DGlobal_600cell.m`，結果：`ans4DGlobal_v4.mx`

| 多胞体 | (C1,C2)総数 | True | False | Unique | 分類 | v3比 |
|--------|:-----------:|-----:|------:|-------:|------|:------:|
| 5胞体  | 20          | 20   | 0     | 20     | **Perfect**   | 同 |
| 8胞体  | 48          | 48   | 0     | 48     | **Perfect**   | 同（geo-unique: 5→1）|
| 16胞体 | 64          | 20   | 44    | 20     | Possible      | 25→20 |
| 24胞体 | 192         | 192  | 0     | 192    | **Perfect**   | 同 |
| 120胞体 | 1,440      | 1,440| 0     | 1,440  | **Perfect**   | **大改善** (1250→1440) |
| 600胞体 | 2,400      | 0    | 2,400 | 0      | **Impossible** | 同 |

**注記**：
- **修正内容**：k≥3 の w タイブレーク時に Det が全候補で ≈0 の場合，xy クロス積 `c2_x·cj_y − c2_y·cj_x` にフォールバック（k=2 と同じ式）。以前はリスト順（非幾何学的）に退化していた
- **8胞体**：セル重心がすべて ±e_i で，k≥3 でも特定の組み合わせで Det=0 が生じリスト順依存 → 修正により geo-unique が 5→1 に（全48 order が幾何学的に等価）
- **120胞体**：Det=0 縮退が頻発していたため 190 ペアが失敗 → 修正により全 1440 ペアが成功し Perfect に
- **16胞体**：25→20 に減少。Det=0 フォールバックがより正確な幾何経路を選ぶため，以前の 25 はリスト順依存の偶発的成功を含んでいたと考えられる
- **600胞体**：修正後も Impossible のまま。正二十面体対称に由来する構造的ボトルネック（全2,400ペアが固定5パターンのステップで停止）であり，タイブレーク修正では解消不可
- k=2 タイブレーク：xy 平面符号付き面積 `c2_x·cj_y − c2_y·cj_x`（+w 軸から見た c2→cj 反時計回り角度）

### 1段バックトラック（`peeling4Df4_1back.m` + RZ，2026-05-11）

計算スクリプト：`run4DBacktrack.m`，結果：`ans4DGlobal1back.mx`

| 多胞体 | (C1,C2)総数 | 1back True | greedy True | 差 | 1back 分類 |
|--------|:-----------:|:----------:|:-----------:|:---:|------|
| 5胞体  | 20          | 20         | 20          | 0   | **Perfect** |
| 8胞体  | 48          | 48         | 48          | 0   | **Perfect** |
| 16胞体 | 64          | 25         | 25          | 0   | Possible |
| 24胞体 | 192         | 192        | 192         | 0   | **Perfect** |
| 120胞体 | 1,440      | 1,250      | 1,250       | 0   | Possible |
| 600胞体 | 2,400      | 0          | 0           | 0   | **Impossible** |

**主な知見**：
- **1段バックトラックはグリーディ版と完全に同一**—どの多胞体でも改善なし
- 600胞体は Impossible のまま（構造的ボトルネックはローカルバックトラックで解消不可）
- アルゴリズム：Det グローバル参照 + RZ 規則（`peeling4Df4.m` 準拠）；justBT フラグで連続バックトラック禁止
- 計算時間：600胞体 516 秒（8.6 分）

### valid net 比率（STLs_v4/，2026-05-12）

計算スクリプト：`run4DExportSTLv4.m`（`ans4DGlobal_v4.mx` から order 読み込み + SAT 重なり判定）

| 多胞体 | Unique orders | Valid (重なりなし) | Overlap | Valid% |
|--------|:-------------:|:-----------------:|:-------:|:------:|
| 5胞体  | 20            | 20                | 0       | 100.0% |
| 8胞体  | 48            | 48                | 0       | 100.0% |
| 16胞体 | 20            | 20                | 0       | 100.0% |
| 24胞体 | 192           | 192               | 0       | 100.0% |
| 120胞体 | 1,440        | 0                 | 1,440   | 0.0%   |
| **合計** | **1,720**  | **280**           | **1,440** | **16.3%** |

**主な知見**：
- 5・8・16・24胞体：全 Apple-Peel order が重なりなし（100%）→ 280枚すべて 3D印刷可能な valid net
- 120胞体：全1,440 order が重なりあり（0%）→ Perfect に改善されても3D空間での valid net は生成しない（構造的問題）
- 600胞体：展開成功が0件（Impossible）のためスキップ
- 計算時間：120胞体 1,440 order × ~0.62 s/order ≈ 896 秒（約15分）
- 旧版（STLs_v3/，285枚）と比べ 16胞体が 25→20 で5枚減，120胞体は Unique 数増加も 0% のまま

### フォールバックなし比較（`peeling4Df4.m` noFB オプション，2026-05-12）

計算スクリプト：`run4DGlobal_nofallback.m`，結果：`ans4DGlobal_nofb.mx`

| 多胞体 | RZ(w) | RZ(n) | RS(w) | RS(n) |
|--------|------:|------:|------:|------:|
| 5胞体  | Perfect/20  | Perfect/20  | Perfect/20  | Perfect/20  |
| 8胞体  | Perfect/48  | Perfect/48  | Perfect/48  | Perfect/48  |
| 16胞体 | Possible/20 | Possible/20 | Impossible/0| Impossible/0|
| 24胞体 | Perfect/192 | Perfect/192 | Perfect/192 | Perfect/192 |
| 120胞体| **Perfect/1440** | **Possible/978** | Impossible/0| Impossible/0|
| 600胞体| Impossible/0| Impossible/0| Impossible/0| Impossible/0|

**主な知見**：
- **フォールバックが有効な唯一の組み合わせは 120胞体 × RZ**：RZ(w)=1440 → RZ(n)=978（462ペア=32%がフォールバック経由）
- 5・8・24胞体では fallback あり/なしで差なし（フォールバックが一切発生しない）
- 16胞体：RZ では fallback の有無で変化なし；RS はフォールバックに関わらず Impossible
- 3D の SnubCube と同じ立ち位置：最も複雑な多胞体でのみフォールバックが機能する
- 600胞体 RZ(n)=275s < RZ(w)=433s：fallback なしで早期終了するため高速（フォールバックが発生していた証拠）
- **フォールバックの幾何学的意味**：120胞体南極付近でセル重心が c1 と w 軸方向に線形従属 → Det≈0 → 全隣接候補が左条件不通過 → fallback（min-w）が正しく「南方向への継続」を提供（3D の真南面選択と同機構）

### RS vs RZ 比較（`peeling4Df4.m`，2026-05-10）

計算スクリプト：`run4DGlobal_RSRZ.m`

| 多胞体 | (C1,C2)総数 | RZ True | RZ 分類 | RS True | RS 分類 |
|--------|:-----------:|--------:|---------|--------:|---------|
| 5胞体  | 20          | 20      | **Perfect**   | 20      | **Perfect**   |
| 8胞体  | 48          | 48      | **Perfect**   | 48      | **Perfect**   |
| 16胞体 | 64          | 25      | Possible      | 0       | **Impossible** |
| 24胞体 | 192         | 192     | **Perfect**   | 192     | **Perfect**   |
| 120胞体 | 1,440      | 1,250   | Possible      | 0       | **Impossible** |
| 600胞体 | 2,400      | 0       | **Impossible** | 0      | **Impossible** |

**主な知見**：
- 4次元では RZ が RS を大幅に上回る
- RS は 16胞体・120胞体で完全失敗（Impossible）
- 5胞体・8胞体・24胞体では RZ と RS が一致
- 3D では RS と RZ が正多面体で完全一致するが，4D では大きな差異が生じる

### 600胞体の終了ステップ分布（RZ，2026-05-10）

計算スクリプト：`analyze_600cell_stuck.m`（全2,400ペアを直接計算）

| 終了ステップ | ペア数 | 割合 |
|:-----------:|-------:|-----:|
| 146         | 742    | 30.9% |
| 150         | 825    | 34.4% |
| 276         | 404    | 16.8% |
| 279         | 323    | 13.5% |
| 284         | 106    |  4.4% |

- 全2,400ペアがわずか5通りの終了ステップのいずれかで停止
- 正二十面体対称に由来する構造的ボトルネックの強い証拠

### 600胞体の cell-centroid-up w 分布と停止ステップ（2026-05-12）

計算スクリプト：`analyze_600cell_wdist.m`，C1=1 で face-up 後の全 600 セル重心 w 座標を解析

#### cell-centroid-up の層構造（31 バンド）

cell-centroid-up では 31 バンドが得られる（Kaino vertex-first の 13 層とは構造が根本的に異なる）．

上から数えた主要な累積境界：

| 層（上から） | バンド数累積 | w 値 | 備考 |
|:------------:|:-----------:|:-----:|------|
| band 9 終端  | 137 | +0.926 | |
| band 10 終端 | 161 | +0.749 | |
| band 15 終端 | 273 | +0.177 | |
| band 16（赤道帯） | 327 | ≈0 | **54 セル（最大バンド）** |

C2 候補（4個）の w 座標は完全一致（spread < 10^{-15}）—600胞体でも同じ定理が成立．

#### 停止ステップと層境界の対応

| 停止ステップ | 最近接（cell-centroid-up） | 差 | 最近接（Kaino） | 差 |
|:-----------:|:-------------------------:|:--:|:---------------:|:--:|
| 146 | 137（band 9 終端） | 9 | 130（layer 4） | 16 |
| 150 | 161（band 10 終端） | 11 | 130（layer 4） | 20 |
| 276 | 273（band 15 終端） | **3** | 270（layer 6） | 6 |
| 279 | 273（band 15 終端） | 6 | 270（layer 6） | 9 |
| 284 | 273（band 15 終端） | 11 | 270（layer 6） | 14 |

**主要な知見**：
- Kaino の vertex-first 層境界との対応という当初仮説は**不支持**（差が 2 倍程度大きい）
- cell-centroid-up 自身の層境界との対応の方が有意に良い
- **ステップ 276-284 ≈ 累積 273**：w≈0.177 の band 15 末尾．直後に**赤道帯（w≈0，54 セル）**が来る → RZ（max-w）で北半球を使い切った後，赤道帯への遷移で詰まる構造的バリアの可能性
- **ステップ 146-150 ≈ 累積 137-161**：w≈0.926〜0.749 の中間帯での別の構造的バリア
- **接続数の制約**：600胞体は隣接数 4（120胞体は 12）→ 各ステップの選択肢が極めて少なく，Det 条件がさらに絞ると行き詰まりやすい

## 3D正多面体の結果（Face-up あり，fallback あり/なし）

| 多面体 | ペア数 | R1(w) | R1(n) | R2(w) | R2(n) | R3(w) | R3(n) |
|--------|:------:|------:|------:|------:|------:|------:|------:|
| Tetrahedron | 12 | 100% | 100% | 100% | 100% | 100% | 100% |
| Cube | 24 | 100% | 100% | 100% | 0% | 100% | 100% |
| Octahedron | 24 | 100% | 100% | 100% | 0% | 100% | 100% |
| Dodecahedron | 60 | 100% | 100% | 100% | 0% | 100% | 100% |
| Icosahedron | 60 | 100% | 100% | 0% | 0% | 100% | 100% |

w = with fallback，n = no fallback
アルゴリズム：大域的 c_1 参照・eps 閾値・単一候補優先なし（2026-05-06 版）
全結果が 0% か 100%—面推移性による等変性理論と一致。

## 3D正多面体の結果（厳密条件 eps < 0, thresh = +10^-10）

計算スクリプト：`run_platonic_strict.m`，結果：`dataPlatonic_strict.mx`（2026-05-07）
RS（Spiral rule）と RZ（Zonal rule）のみを計算。

| 多面体 | ペア数 | RS(w) | RS(n) | RZ(w) | RZ(n) |
|--------|:------:|------:|------:|------:|------:|
| Tetrahedron | 12 | 100% | 100% | 100% | 100% |
| Cube | 24 | 100% | 0% | 100% | 0% |
| Octahedron | 24 | 100% | 0% | 100% | 0% |
| Dodecahedron | 60 | 100% | 0% | 100% | 0% |
| Icosahedron | 60 | 100% | 0% | 100% | 0% |

- **with fallback は常に 100%**（「真南面」が fallback 経由で選ばれる）
- **no fallback は Tetrahedron のみ 100%**（正四面体には「真南面」が存在しない）
- **RS と RZ が完全一致**（正多面体の緯度帯構造による）

## 3Dアルキメデス多面体の結果（Face-up あり，fallback あり/なし）

計算スクリプト：`run_archimedean_faceup.m`（fallback あり），`run_archimedean_nofallback.m`（fallback なし，2026-05-12再実行）
結果：`archimedean_faceup_results.mx`，`archimedean_nofallback_results.mx`
アルゴリズム：大域的 c_1 参照・eps 閾値・単一候補優先なし・**Det ベース選択**（2026-05-12版）
RS 選択：`leftCands` から min Det（最大右旋回），RS フォールバック：Darboux frame min φ
RZ 選択：`leftCands` から max z → min Det タイブレーク，RZ フォールバック：min z → min Det

| 多面体 | F | Pairs | RS(w) | RS(n) | RZ(w) | RZ(n) |
|--------|:-:|------:|------:|------:|------:|------:|
| TruncatedTetrahedron | 8 | 36 | 24 | 24 | 12 | 12 |
| Cuboctahedron | 14 | 48 | 0 | 0 | 0 | 0 |
| TruncatedCube | 14 | 72 | 0 | 0 | 0 | 0 |
| TruncatedOctahedron | 14 | 72 | 30 | 30 | 72 | 72 |
| Rhombicuboctahedron | 26 | 96 | 0 | 0 | 0 | 0 |
| TruncatedCuboctahedron | 26 | 144 | 0 | 0 | 144 | 144 |
| **SnubCube** | 38 | 120 | **24** | **0** | **48** | **24** |
| Icosidodecahedron | 32 | 120 | 0 | 0 | 0 | 0 |
| TruncatedDodecahedron | 32 | 180 | 0 | 0 | 0 | 0 |
| TruncatedIcosahedron | 32 | 180 | 0 | 0 | 180 | 180 |
| Rhombicosidodecahedron | 62 | 240 | 0 | 0 | 0 | 0 |
| TruncatedIcosidodecahedron | 62 | 360 | 0 | 0 | 240 | 240 |
| SnubDodecahedron | 92 | 300 | 0 | 0 | 0 | 0 |

w = with fallback，n = no fallback

主な知見：
- **7種（54%）が RS・RZ とも 0%**（Cuboctahedron, TruncatedCube, Rhombicuboctahedron, Icosidodecahedron, TruncatedDodecahedron, Rhombicosidodecahedron, SnubDodecahedron）
- **RZ > RS の4種**：TruncatedOctahedron（RZ=100% > RS=41.7%），TruncatedCuboctahedron（RZ=100% > RS=0%），TruncatedIcosahedron（RZ=100% > RS=0%），TruncatedIcosidodecahedron（RZ=66.7% > RS=0%）
- **RS > RZ の1種**：TruncatedTetrahedron（RS=66.7% > RZ=33.3%）
- **RS = RZ は 8/13種（いずれも 0%）**；RZ は 4 種で 100% を達成するが RS はゼロ以外の上限は 66.7%
- **RZ 結果は Det ベース選択導入前後で完全不変**（min-φ タイブレーク ≡ min-Det タイブレーク）
- **RS 結果は大幅に変化**（min-Det ≠ min-φ）：旧版では TruncatedCuboctahedron/TruncatedIcosahedron/TruncatedIcosidodecahedron も RS=100% だったが，Det ベースでは 0%
- **フォールバック効果（2026-05-12）**：RS・RZ とも 12/13 種でフォールバック除去による変化なし．SnubCube のみ例外（RS: 24→0，RZ: 48→24）
- **等変性**：各成功数は面軌道ごとに一様（RZ は旧版と変わらず13種すべて等変；RS は面タイプが2種以上ある多面体では部分的成功（例 TruncatedOctahedron 30/72=41.7%）も等変と矛盾しない）
- **鏡像対称性（2026-05-10 確認済み，RZ のみ）**：全13種で RZ の成功数が鏡像と一致；キラルの SnubCube は成功ペアが鏡像間で入れ替わる

---

## 120胞体 RZ nets の幾何学的多様性（C1=3，2026-05-14）

スクリプト：`analyze_c2_xyz.m`（セッション内で作成・実行）

### Geo-unique net 数と軌道構造

C1=3，cell-centroid-up で 12 通りの C2 すべてが成功するが，`netGeoKey4D` による重複除去で **7 種類の幾何学的に異なる net** しか存在しない。

| 軌道 | C2 メンバー | サイズ |
|------|:----------:|:------:|
| orbit 1 | 1, 15, 17, 97 | 4 |
| orbit 2 | 4, 18, 51 | 3 |
| singleton × 5 | 2, 10, 19, 100, 115 | 各 1 |

確認方法：`ans4DGlobal_v4.mx` から C1=3 の全 order を取得し，geo-key でグルーピング。

### C2 重心の正二十面体配置

Cell-centroid-up 後の 12 個の C2 セル重心（xyz 部分）は，黄金比 φ≈1.618 を使った循環座標
$(0, \pm\varphi/2, \pm\varphi^2/2)$ の全置換により正十二面体の面配置に対応した**正二十面体**を形成する
（辺長 φ，各 r_xyz≈1.539）。これは C1 が正十二面体セルであることの直接的帰結。

### 対蹠ペア

12 個の C2 は6対の対蹠ペアを形成：{1,18},{2,17},{4,19},{10,115},{15,100},{51,97}。
対蹠ペアは k=2 タイブレーカー（xy クロス積スコア）の符号が逆転するため，**同じ軌道に属さない**。

### 普遍的内部コア（k=1..13）

7 種の net はすべて，最初の 13 セル（C1 + 12 個の直接隣接セル）について
r_3D（3D 展開図内の重心距離）が完全に一致する：

```
r_3D[k=1..13] = {0.00, 1.70, 1.79, 1.80, 1.81, 1.82, 1.81, 1.96, 1.99, 2.04, 2.06, 2.02, 2.09}
```

k=14 以降で分岐が始まる。

### 螺旋構造パターン（累積旋回角による分類）

各 net の円柱座標 (r_xy, z, θ) を展開順に計算し，累積旋回角（単位：周回数）を比較：

| C2 | 軌道 | 旋回数 | 方向 | z 範囲 | パターン型 |
|----|:----:|:------:|:----:|:-------:|:------:|
|  1 | {1,15,17,97} | +1.07 | CCW | −3〜+26 | A（CCW 上昇） |
| 10 | singleton | +0.85 | CCW | −5〜+8  | A（CCW 上昇） |
| 19 | singleton | +1.08 | CCW | −4〜+20 | A（CCW 上昇） |
|100 | singleton | +0.75 | CCW | −3〜+29 | A（CCW 上昇） |
|115 | singleton | +0.82 | CCW | −13〜+6 | A（CCW 下降） |
|  2 | singleton | −1.14 | CW  | −7〜+4  | B（CW 反転） |
|  4 | {4,18,51} | −0.07 | ≈0  | −4〜+23 | C（柱状） |

- **型 A**（5/7 net）：RZ の max-xy-cross-product タイブレーカーが CCW を選好することと整合
- **型 B**（1/7，C2=2）：CW の反転型。C2=2 は C2=17 と対蹠（C2=17 は orbit 1 で CCW）
- **型 C**（1/7，C2=4 代表）：旋回ほぼゼロ，z 方向に伸びる柱状構造。orbit {4,18,51} に対応

### 注意：`Round[x, N]` バグ

Mathematica の `Round[x, N]` は N 小数点以下ではなく **N の最近傍倍数** に丸める。
3 桁表示には必ず `Round[x, 0.001]` を使うこと（`Round[x, 3]` は ×3 の倍数に丸まる）。

### 120胞体の自己交差構造解析（2026-05-14）

スクリプト：`analyze_120cell_overlap_bands.m`，`analyze_120cell_band_backtrack.m`

#### 帯間 vs 帯内の重なり分布

C1=3 の全 12 C2 × 各 308 重なりペアを解析した結果：
- **同一帯内の重なり：88 ペア（28.6%）**
- **異なる帯間の重なり：220 ペア（71.4%）**

異帯重なりはほぼすべて **隣接帯間（帯k ↔ 帯k+1）** であり，遠く離れた帯の衝突ではない。
自己交差は帯境界で各セルが隣接帯に「はみ出す」ことで発生している。

#### 帯境界代替エントリー検証（2026-05-14）

仮説：帯遷移ステップで RZ が選ぶセル以外を選べば重なりを回避できるか？

スクリプト `analyze_120cell_band_backtrack.m` を C1=3，C2=1 で実行した結果：
- **全帯遷移ステップで「代替なし」（隣接する次帯セルが1個のみ）**
- RZ アルゴリズムは帯境界を越える際に選択の余地がない

追加の観察：
- RZ の帯走査は単純な「帯を順番に一周」ではなく，帯4↔5，5↔6，6↔7，7↔8 の間で何度もジグザグする
- 帯遷移が約 60 回（全 119 ステップ中）発生し，各遷移で唯一の通路しか存在しない

**結論**：帯境界での局所的なセル置換で重なりは回避できない。重なりを解消するには，帯内の巡回順序（どのセルを最後に置くか）を変えるグローバルな最適化が必要。これはグリーディ系では解決困難。

---

## Face-up 方式の比較実験（120胞体 C1=3，2026-05-11）

スクリプト：`analyze_120cell_c2types.m`，`test_faceup_pair.m`，`test_faceup_allC2.m`

### 幾何学的発見：C2 候補の w 座標が全て同値

120胞体の任意の C1 について，face-up 後（cell-centroid-up）の C2 候補 12 個の w 座標の
spread は ~10^{-15}（浮動小数点の機械イプシロン）で，**実質的に完全一致**。

```
C1=3, 10, 50, 60 で確認：w-spread < 10^{-15}（全12候補）
```

**数学的理由**：C1 の安定化部分群（≅ I_h，位数120）は 12 個の C2 候補を推移的に置換し，
+w 軸（face-up 方向 = C1 重心方向）を保存する。よって C2 全員が同じ w 投影値を持つ。
これは 120胞体に限らず**任意の正多胞体の任意の C1** で成立する定理。

**含意**：RZ 規則の "max w" 基準は k=2 では常にタイ → C2 の選択は k=2 タイブレーカー
（xy 平面反時計回り角度スコア）のみで決まる。

### Face-up 方式の比較

| 方式 | C2 w-bands | spread |
|------|:---------:|-------:|
| cell-centroid-up（現行）| 1（全同値）| ~10^{-15} |
| vertex-up | 4（3+3+3+3）| 0.926 |
| 3D-face-centroid-up | 4（**1+5+5+1**）| 0.951 |

**3D-face-centroid-up**：C1 と C2 の共有ペンタゴン面の重心を +w に揃える。
3D の face-up（開始面の重心を +z）の**直接類推**であり，正十二面体の 1+5+5+1 緯度帯構造が復元される。

### C1=3，全 C2（12通り）での比較結果（v4 fix 適用後，2026-05-13）

**注記**：旧テスト（2026-05-11，10通り）は `peel1Pair` に k≥3 Det=0 フォールバックが未実装で，
C2=19・C2=51 が漏れていた。`gen_120cell_faceup_nets.m` に v4 fix を適用し全 12 C2 を再計算。

| C2 | A: cell-centroid-up | A-ov | B: 3D-face-centroid-up | B-ov | B−A |
|----|:-------------------:|:----:|:----------------------:|:----:|:---:|
|  1 | OK | 27 | OK | **14** | −13 |
|  2 | OK | 25 | OK | 23 | −2 |
|  4 | OK | 21 | FAIL | — | — |
| 10 | OK | 26 | FAIL | — | — |
| 15 | OK | 27 | FAIL | — | — |
| 17 | OK | 27 | FAIL | — | — |
| 18 | OK | 21 | OK | 22 | +1 |
| 19 | OK | 25 | OK | **15** | −10 |
| 51 | OK | 21 | OK | 27 | +6 |
| 97 | OK | 27 | OK | **13** | −14 |
|100 | OK | 26 | OK | **16** | −10 |
|115 | OK | 35 | FAIL | — | — |
| **合計** | 12/12 成功 | **308** | 7/12 成功 | **130** | — |

**知見**：
- **B 成功（7ケース）**：C2 = 1, 2, 18, 19, 51, 97, 100
  - B が A より改善（ov 減少）：C2 = 1(−13), 2(−2), 19(−10), 97(−14), 100(−10)
  - B が A より悪化（ov 増加）：C2 = 18(+1), 51(+6)
- **B 失敗（5ケース）**：C2 = 4, 10, 15, 17, 115
- 旧テストとの主な差分：C2=115 が OK→FAIL，C2=2/19/51/100 が FAIL→OK（いずれも v4 fix による経路変化）
- 成功7ケース中5ケースで重なり削減（平均改善幅 −9.8），2ケースで悪化

### 3D-face-centroid-up を 600胞体に適用した結果（2026-05-11）

スクリプト：`test_600cell_faceup.m`（全 2400 ペア，227 秒）

**結論：600胞体は 3D-face-centroid-up でも Impossible（0/2400）**

| 方式 | 成功数 | 停止ステップ種類 | 主な停止ステップ |
|------|:------:|:--------------:|----------------|
| cell-centroid-up（現行）| 0/2400 | **5通り** | 146, 150, 276, 279, 284 |
| 3D-face-centroid-up     | 0/2400 | **25通り** | 88〜201 に分散 |

- 停止ステップが 5通り → 25通りに多様化（対称性ボトルネックのパターンが変化）
- 最大停止ステップが 284 → 201 に短縮（より早く詰まる）
- 主な停止点: step 99（320件, 13.3%），step 103（319件, 13.3%）

**解釈**：対称性は破れてボトルネックのパターンは変化したが，600胞体の Impossible は不変。
停止ステップの多様化は「1つの深い障壁」→「複数の浅い障壁」への変化を示唆。
RZ 規則が 600胞体の構造に根本的に合っていない可能性が高い。

### vertex-up を 600胞体に適用した結果（2026-05-12）

スクリプト：`test_600cell_vertexup.m`（C1=1，4頂点 × 4 C2 候補 = 16 ケース）

**結論：vertex-up でも全 16 ケースが STUCK（成功 0）**

| 方式 | 成功数 | 停止ステップの範囲 | 備考 |
|------|:------:|:-----------------:|------|
| cell-centroid-up（現行）| 0/4 | 146, 284 | 2パターン（正二十面体対称） |
| vertex-up（v1〜v4 × C2） | 0/16 | 80〜206 | 14種類に分散，多くが早く詰まる |

- 停止ステップは 80〜206 に散らばった（cell-centroid-up の 5 固定パターンとは完全に異なる）
- 多くのケースで cell-centroid-up より**早く**詰まる → vertex-up は「方向の多様性」をもたらすが，w-based 道案内が不正確になり逆効果
- face-centroid-up と同じパターン：「詰まる場所は変わるが，詰まること自体は変わらない」

**グリーディ系 up 方向変更の網羅的まとめ**：

| up 方向 | 規模 | 成功数 | 停止パターン |
|---------|:----:|:------:|------------|
| cell-centroid-up（現行） | 2400 ペア | 0 | 5 固定（正二十面体対称） |
| 3D-face-centroid-up | 2400 ペア | 0 | 25 種類に分散 |
| vertex-up | 16 ケース（C1=1） | 0 | 14 種類に分散，平均的に早く詰まる |

**結論**：up 方向の変更は投影構造を変えるが，「接続数 4 という貧しいグラフ構造の中で Det 条件が全候補を封鎖する」という根本的な障壁を解消できない。
600胞体の Impossible はグリーディ系アルゴリズムの構造的限界であり，up 方向・規則の組み合わせを変えても打破できない可能性が高い。

### 未解決の問い（次のセッションへの引継ぎ）

**Q1（v4 fix 後更新）: 120胞体で B が成功する C2（1, 2, 18, 19, 51, 97, 100）と失敗する C2（4, 10, 15, 17, 115）に幾何学的な共通点はあるか？**

旧仮説（C2=115 が OK だった時代の4バンド仮説）は v4 fix で成功セットが変わったため再検討が必要。
旧バンド仮説は無効（C2=115 は最高バンドだったが FAIL，C2=2/19/51/100 は中位バンドで OK）。

**Q2（検討終了）：600胞体の Impossible を打破する up 方向・規則の組み合わせはあるか？**
- cell-centroid-up + RZ：0/2400
- 3D-face-centroid-up + RZ：0/2400（停止パターンのみ変化）
- vertex-up + RZ：0/16（早く詰まる）
- RS 規則：0/2400（cell-centroid-up で確認済み）
- **結論**：グリーディ系では打破不可能と判断。多段バックトラック・探索的アプローチが必要だが，計算コストが膨大で現実的でない。

---

## 今後の課題（優先度順）

| 優先度 | タスク |
|:------:|--------|
| 済 | RS（螺旋規則）を `peeling4Df4.m` に実装し，全6胞体で RZ と比較（`run4DGlobal_RSRZ.m`，2026-05-10）：4Dでは RZ が RS を大きく上回り，RS は 16胞体・120胞体で Impossible |
| 済 | 600胞体の終了ステップ頻度分布を解析（`analyze_600cell_stuck.m`，2026-05-10）：全2,400ペアがわずか5通りのステップ（146,150,276,279,284）で停止 |
| 済 | 1段バックトラック（`peeling4Df4_1back.m`，Det グローバル参照 + RZ）で全6胞体計算（`run4DBacktrack.m`，2026-05-11）：**グリーディ版と完全に同一**—改善なし（600胞体は Impossible のまま，120胞体 1250/1440，16胞体 25/64） |
| 済 | k≥3 Det=0 縮退のとき xy クロス積フォールバックを実装（`peeling4Df4.m` 修正，2026-05-12）：8胞体 geo-unique 5→1，**120胞体 1250→1440（Perfect）**，600胞体は Impossible のまま |
| 済 | 4D正多胞体のフォールバックなし計算（`peeling4Df4.m` noFB オプション追加，`run4DGlobal_nofallback.m`，2026-05-12）：**120胞体 RZ のみ** fallback が有効（1440→978），他は変化なし |
| 済 | STL valid net 再集計（`run4DExportSTLv4.m`，2026-05-12）：5/8/16/24胞体 100%，120胞体 0%，合計 280/1720（16.3%）→ `STLs_v4/` に保存 |
| 済 | k=2 タイブレークの幾何学的基準を実装（2026-05-10）：xy 平面符号付き面積 c2_x·cj_y - c2_y·cj_x（+w 軸から見た c2→cj 反時計回り角度）。k≥3 は max Det，k=2 は xy クロス積で一貫した幾何学的タイブレーク |
| 済 | xy クロス積 k=2 タイブレーク版で全C1計算し直す（`run4DGlobal_xy2.m` 実行 → `ans4DGlobal_v3.mx`，2026-05-10）：16胞体のみ変化（34→25），他は同値 |
| 済 | `peeling4Df4.m` を新規作成：c1–c2 平面グローバル参照 `Det[{c1,c2,c_k,c_j}] >= -eps`（2026-05-08） |
| 済 | `run4DPeelingGlobal.m` で全6胞体を計算・比較（2026-05-08）：5/8胞体 Perfect，24/120胞体改善，16胞体悪化，600胞体 Impossible |
| 済 | `peeling4Df4.m` に max-w + Det タイブレーク実装，`run4DGlobal_tiebreak.m` で全C1計算（2026-05-10）：16胞体 34/64，24胞体 Perfect，120胞体 1250/1440 |
| 済 | アルキメデス多面体のフォールバックなし計算：`run_archimedean_nofallback.m` を実行，13種中12種で変化なし，SnubCube のみ 40.0%→20.0%（2026-05-10） |
| 済 | 3D-face-centroid-up を 600胞体に適用（`test_600cell_faceup.m`，2026-05-11）：Impossible のまま（0/2400），停止ステップが 5種類 → 25種類に多様化，最大 284→201 に短縮 |
| 済 | `gen_120cell_faceup_nets.m` に v4 k≥3 Det=0 fix 適用・全12 C2 対応・STL 出力（`STLs_120cell_faceup/`，2026-05-13）：face-up 成功 C2 が {1,18,97,115} → {1,2,18,19,51,97,100}（7ケース）に更新 |
| 済 | `paper_draft.tex` の face-up 比較テーブル・本文を全面更新（12 C2，7成功，avg 24.6→18.6，2026-05-13） |
| 済 | `paper_draft.tex` 全体レビュー・投稿準備（2026-05-14 後半）：下記「論文草稿の現状」参照 |
| 済 | `unfold4D.m` grayscale 色（0.2→0.7）・黒エッジ・Neutral lighting・コメント修正（2026-05-13） |
| 済 | `gen_snubcube_both.m` 変数名バグ修正（アンダースコア→camelCase）・Rasterize 出力（2026-05-13） |
| 済 | `visualize_archimedean_nets.m` に SnubCube mirror（右手系）inline 表示を追加（2026-05-13） |
| 済 | `paper_draft.tex` Figure 6 を `snubcube_both.pdf` に更新，Figure 7 重複削除，4D nets figure 追加（2026-05-13） |
| 済 | vertex-up を 600胞体に適用（`test_600cell_vertexup.m`，2026-05-12）：全16ケース STUCK．**グリーディ系では up 方向変更で打破不可と結論** |
| 済 | 4D正多胞体のフォールバックなし計算：`run4DGlobal_nofallback.m`（2026-05-12）：120胞体 RZ のみ fallback が有効（1440→978），他は変化なし |
| 済 | STL valid net 比率集計（`run4DExportSTLv3.m`，2026-05-11）→ `STLs_v3/`（285枚）；v4 版（`run4DExportSTLv4.m`，2026-05-12）→ `STLs_v4/`（280枚） |
| 済 | `peeling4Df4_1back.m` グリーディと完全同値（2026-05-11） |
| 済 | `peeling3DLoxo.m` 大域 c_1 参照・eps・単一候補優先廃止（2026-05-06）；厳密条件計算（2026-05-07）；アルキメデス再計算（2026-05-07）；SnubCube バグ修正（2026-05-08） |
| 済 | 600胞体 w 分布解析（`analyze_600cell_wdist.m`，2026-05-12）：31 バンド，停止ステップは cell-centroid-up 層境界に対応 |
| 検討済・却下 | Face-up 方式の全面切り替え（2026-05-12）：600胞体は Impossible のまま，他は改善なし → 実施しない |

### 未解決・未実施タスク

| 優先度 | タスク |
|:------:|--------|
| 低 | Kaino (2019) との内容面の比較（未実施） |

### 検討済み・対象外

| タスク | 理由 |
|--------|------|
| Q1（face-up 成功/失敗 C2 の幾何学的違い） | face-up 削除と同じ理由でスコープ外（2026-05-14） |
| orbit {4,18,51} の柱状構造の幾何学的理由 | 主結論を変えない（2026-05-14） |
| 600胞体の構造的ボトルネックへの深い対処 | 構造的理由は Discussion に記載済み；多段 BT は非現実的（2026-05-14） |
| 2段階選択（RZ→R2）の実装と評価 | 論文の見通しを悪くする（face-up 削除と同じ理由） |
| 左利きペーリングとの対称性検討 | 結果が鏡像になるだけの可能性が高い |
| 4D nets figure（`4d_nets.pdf`）の生成 | 済（2026-05-14） |
| face-up セクションの最終判断 | 済・削除（Discussion に要約として移動，2026-05-14） |

---

## 論文草稿の現状（`paper_draft.tex`，2026-05-16 更新）

### 投稿準備状況

- **arXiv 投稿は即可**
- **CGTA（第1候補）は submittable，ただし major revision を覚悟すべき水準**
- **DCG レベルを狙うなら追加の形式化（命題化・600胞体の structural argument）が必要**

詳細は下記「Opus 4.7 レビュー（2026-05-16）」および `commentsOpus260516.md` を参照。

### Opus 4.7 レビュー（2026-05-16）

#### 投稿水準に達している点

| 項目 | 評価 |
|------|------|
| 新規性の明示 | ✓ Kaino, Akitaya との差別化が明確（4D 統一グリーディ + 全6胞体・全 starting pair） |
| アルゴリズム提示 | ✓ Algorithm 1, 2 で 3D/4D とも疑似コード完備 |
| 実験スコープ | ✓ 24 立体（5+13+6）×2 規則 × fallback 有無 |
| 等変性の議論 | ✓ Companion との対比で動機付け明快 |
| ordering vs validity の区別 | ✓ Abstract レベルで 120-cell を honest に説明 |
| Mirror symmetry / chirality | ✓ 副次的結果として位置付け良好 |
| 限界の議論 | ✓ 600-cell を multiple variants で検証し Impossible を確認 |
| 図表・参照整合性 | ✓ 直近のコンパイルで警告なし、未定義参照なし |

#### 査読で指摘される可能性が高い点（major revision 級）

1. **等変性結果が「Remark」止まり**：Section 2 の等変性議論は数式は揃っているが正式に命題化されていない。CGTA・DCG では `Proposition` + `Proof` 形式にすべき。
2. **600-cell の "icosahedral bottleneck" が経験則止まり**：5 通りの停止ステップ固定（146, 150, 276, 279, 284）の構造的説明が薄い。各停止ステップで Det 条件が全候補を排除する具体例を1つ示すか，セル隣接グラフの構造的補題が必要。
3. **コード/データ可用性ステートメントの欠如**：Mathematica 実装（`peeling3DLoxo.m`, `peeling4Df4.m`, `_4DData/`）の公開予定を明記すべき。近年の computational geometry 系誌は事実上必須。
4. **参考文献が9件と寡少**：Schlickenrieder（spiral unfolding の起源），Pak，Bern et al.（Ununfoldable polyhedra），4D 関連の Coxeter, Towle 等を追加検討。
5. **「なぜこの2規則か」の動機付けが弱い**：他の plausible rule（min angle, nearest neighbor 等）と比べてなぜ RS, RZ を選んだかを1段落追加すべき。

#### minor だが対応すれば質が上がる

- Section 5（Computational Examples）が Section 3, 4 と一部重複 — 個別深掘りに整理
- Conclusion の Future Work が短い — より具体的に
- Abstract に Perfect/Possible/Impossible 分類が貢献の一部であることを明記
- Acknowledgements に ORCID, conflict-of-interest 宣言（誌指定で必要になる場合あり）

#### 推奨アクション

| 投稿先 | 必要な追加作業 | 期間目安 |
|--------|---------------|---------|
| arXiv（即時） | なし | 即日 |
| **CGTA**（第1候補） | 上記 1–3 を対応 | 1–2 週間 |
| DCG（高目標） | 上記すべて + 600-cell の structural argument | 1–2 ヶ月 |
| J. Comput. Geom. | 1, 3 のみ対応で十分 | 数日 |

**推奨フロー**：arXiv に先に置いて DOI を確保し，並行して CGTA 向けに 1, 3 だけ最低限対応してから submit。3 (code availability) は GitHub レポジトリ準備で 1 日，1 (Proposition 化) は既存議論の整形で半日〜1日。

### このセッション（2026-05-16）で行った修正

| 修正内容 | 詳細 |
|----------|------|
| line 671: "four" → "three" | RZ 100% は3種（TruncOcta, TruncCubocta, TruncIcosa） |
| line 671–674 文構造整理 | "(all 100%)" 削除，"and is best on" → 完全な文に |
| line 703: "four solids achieving 100%" → "three solids" | Truncated Icosidodecahedron は 66.7% なので除外，補足文追加 |
| line 1258–1259: "verified computationally for all four polytopes" | 16/24-cell の geo-unique=1 確認に基づく根拠追加 |
| line 1386: "9 symmetric bands" → "9 bands (grouped by symmetry orbit)" | symmetry orbit による粗い帯であることを明示 |
| line 1394: `$w \approx 2.118$` 削除 | 正規化が未定義の数値を削除 |
| Conclusion: SnubCube/TT の partial 結果追加 | 13種の内訳（3+1+2+7=13）が完結するよう修正 |
| 4D 条件式に `\label{eq:4dleft}` 付与 | `\[...\]` → `equation` 環境に変更 |
| line 1094–1095: `eq:left`（3D）→ `eq:4dleft`（4D），"base point" → "reference plane (c₁,c₂)" | 4D 版の正しい説明に修正 |
| Figure 5 キャプション | TT は RZ で表示中だが RS が高い（66.7%）ことを明記 |
| `tab:3dnet` キャプション | 600-cell の em-dash が Impossible を意味することを `tab:4dglobal` 参照付きで補足 |

### 2026-05-14 後半の修正

| 修正内容 | 詳細 |
|----------|------|
| Section 4.2 削除 | "4D Darboux Frame via xy-Projection"（`peeling4Df.m` 旧版の説明，主アルゴリズムでは不使用） |
| Table 4dcomp 削除 | xy-projection vs global reference の比較表（Section 4.2 削除に伴い） |
| Observations 内の xy-projection 言及を削除 | 24胞体・120胞体・16胞体の Observations から旧数値（119/192，357，26）を除去 |
| "three-rule" → "two-rule" | 序文 line 139 |
| Dodecahedron "RS" → "any" | tab:summary（RS も RZ も 100% のため） |
| 5件の未引用文献に引用追加 | Demaine2007（序文），AronovORourke1992（序文），Akitaya2024（序文 Kaino 文に併記），Buekenhout1998・Devadoss2022（4D Setup 節末） |
| `\ref{sec:4dresults-global}` 修正 | 未定義ラベル → `\ref{sec:4dresults}` |
| Dodecahedron 節の誤図参照削除 | fig:peel-example（Truncated Icosahedron 図）を Dodecahedron 節で参照していた文を削除 |
| Abstract に ordering vs validity 区別を明記 | 120胞体が Perfect でも 3D valid net は 0 件であることを Abstract レベルで説明 |
| Perfect/Possible/Impossible 定義を前出し | Section 3 冒頭（以前は Discussion にのみあった） |
| Discussion 新節追加 | "Spiral vs. Zonal Strategy in Four Dimensions"：4D で RS が劣化する構造的理由を説明 |

### 投稿候補ジャーナル

| 優先度 | ジャーナル | 備考 |
|:------:|-----------|------|
| 第1候補 | Computational Geometry: Theory and Applications（Elsevier） | Devadoss2022 の掲載誌，内容適合度最高，Scopus 収録確実 |
| 第2候補 | Discrete & Computational Geometry（Springer） | 分野最高峰，競争率高い |
| 第3候補 | Journal of Computational Geometry | オープンアクセス，Scopus 収録要確認 |
| 代替 | Graphs and Combinatorics（Springer） | 等変性・分類の理論面を前面に出す場合 |

### 論文の現在の構成（2026-05-16 更新）

```
1. Introduction
2. Algorithm
   2.1 Geometric Setup
   2.2 Two Selection Rules
   2.3 Net Construction from Selection Order
3. Results: 3D Polyhedra（Perfect/Possible/Impossible 定義はここで）
   3.1 Platonic Solids
   3.2 Archimedean Solids
   3.3 Mirror Symmetry and Chirality Detection
   3.4 Truncated Icosahedron: Hexagonal-Band Path under RZ
   3.5 Snub Cube: Equivariance and Chiral Structure
4. Extension to Four Dimensions
   4.1 4D Algorithm
   4.2 Results: Regular 4-Polytopes
   4.3 Three-Dimensional Realisation of 4D Unfoldings
5. Computational Examples
   5.1 Cross-Dimensional Summary（tab:summary）
   5.2 5-Cell: Agreement of RS and RZ
   5.3 16-Cell: RZ Stalling
   5.4 120-Cell: Perfect Result under RZ
6. Discussion
   [無題段落] Role of Face-Type Uniformity
   6.1 Spiral vs. Zonal Strategy in Four Dimensions
   [無題段落] The 600-Cell
   6.2 Relation to the Companion Implementation
7. Conclusion
```

### 注意事項

- bibitem キーは `\bibitem{Yoshino2026arXiv}`（arXiv 番号 2604.16204，2026年）
- 3D RS フォールバック（Darboux Frame）は Algorithm 節本文で説明，4D アルゴリズムには使用しない
- Section 4.2 で説明していた xy-projection アプローチ（`peeling4Df.m`）は **論文から削除済み**。比較データも削除済み

---

## 対称性に関する理論的考察

### 正多面体における面推移性とアルゴリズムの等変性

正多面体は**面推移的**（face-transitive）であり，任意の2面を写す回転対称が存在する。
ここで，アルゴリズムが対称群 G の下で**等変**（equivariant）かどうかを考える。

**結論：アルゴリズムは G の下で等変であり，成功率は 0% か 100% のいずれかでなければならない。**

証明の概略：σ ∈ G が (F1, F2) → (F1', F2') を与えるとき，

1. **Face-up 回転後の重心**：`cc_{F1'}[σ(j)] = A · cc_{F1}[j]`
   ただし A = Face-up(F1') ∘ σ ∘ Face-up(F1)⁻¹ は proper rotation（det = 1）
2. **左半空間条件**：`Det[{A·a, A·b, A·c}] = det(A)·Det[{a,b,c}] = Det[{a,b,c}]`
   → 左候補集合は σ で置換されるだけで，条件の成否は変わらない
3. **R1/R2（Det 選択）**：`Det[{A·c_1, A·c_k, A·c_j}] = Det[{c_1,c_k,c_j}]`（det(A)=1）→ min/max Det の順序は σ で不変
4. **R3（max z）**：z 座標は σ によって面ラベルが置換されるだけで，
   同じ面（σ でリラベルされた）が選ばれる

したがって，(F1, F2) で成功 ⟺ (F1', F2') で成功。
面推移的な正多面体では全ペアが等価なので，成功率は 0% か 100% のみ。

### 数値誤差問題の解決（2026-05-06）

旧実装では Dodecahedron R3(with fallback) = 53.3% など，等変性理論（0% or 100%）に矛盾する結果が生じていた。原因と対処：

**原因**：
1. 左条件の参照点が局所的（c_{k-1}）だったため，螺旋方向が各ステップで変化し対称性が崩れていた
2. 正十二面体の黄金比座標を `N[...]` で有限精度化すると，理論上0となるはずの行列式値が ~10⁻¹⁶ の誤差を持ち，符号判定がペアによって逆転していた
3. 単一候補優先が左条件をバイパスし，大域螺旋保証を破っていた

**対処**：
- 左条件を `Det[{cc[[top]], cc[[last]], cc[[j]]}] >= -eps`（c_1 大域固定参照・`eps = 10^-10`）に変更
- 単一候補優先を廃止（全候補数に関わらず左条件を適用）

**結果**：全5種の正多面体で全結果が 0% か 100% になり，等変性理論と整合。

---

## コンパニオン論文（arXiv:2604.16204）との実装差異

`peeling3DLoxo.m`（本論文）はコンパニオン論文の `peeling3Df` と以下の3点で異なる。

| 項目 | 旧コード `peeling3Df` | 現行 `peeling3DLoxo.m` |
|------|----------------------|------------------------|
| フィルタ参照点 | c_k（現在の面，毎ステップ変化），厳密 `> 0`，ε なし | c_1（開始面，固定），`<= ε`，ε = 10^{-10} |
| タイブレーク | リスト順（`Position` の最初の要素）— 非幾何的 | min Det — 幾何学的・等変 |
| RS 規則 | なし | あり（min Det 選択） |

旧コードの c_k 参照とリスト順タイブレークは**暗黙的に等変性を破る**。
フィルタ方向が毎ステップ回転し，タイ結果が面の番号付け順（`PolyhedronData` 内部順）に依存するため，
回転対称操作で (F1, F2) が別ペアに写されると同じ結果が保証されない。
結果として正十二面体で 53.3%（等変性理論は 0% か 100% のみを予測）という矛盾した値が現れていた。

現行実装の c_1 固定参照・min Det タイブレーク・ε 閾値が等変性を回復し，数値的にも堅牢。

---

## コードを読むときの注意

- `summary.tex`（または `summary_en.tex`）が最も包括的な情報源。コードの意図が不明な場合はまずここを参照。
- `peeling4Df.m` と `peeling4Df3.m` の**左条件の実装が異なる**点に注意（xy 2成分 vs xyz 3成分行列式）。
- 3D版の右条件（フィルタ）は `Det[{cc[[top]], cc[[last]], cc[[j]]}] <= eps`（c_1 大域固定参照，`eps = 10^-10`）。`peeling3DLoxo.m` の `p3lPeelPair` を参照。
- 3D版の**選択基準も Det ベース**（2026-05-12）：RS は min Det，RZ は max z → min Det タイブレーク。RS フォールバックのみ Darboux frame min φ。
- **単一候補優先は廃止**（2026-05-06）：大域螺旋保証のため，単一候補でも右条件を適用する。4D版（`peeling4Df.m`）は旧実装のままなので注意。
- アルキメデス多面体の計算結果（`archimedean_faceup_results.mx`, `archimedean_nofallback_results.mx`）は 2026-05-12 に Det ベース選択で再計算済み。

---

## 可視化

### `visualize_platonic_nets_v2.m`（現行版）

```mathematica
Get["/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/visualize_platonic_nets_v2.m"]
```

`dataPlatonic.mx`（標準）と `dataPlatonic_strict.mx`（厳密）を読み込み，
5種の正多面体 × RS/RZ × withFallback/noFallback の展開図を描画する．
R2 (loxodrome) は廃止済みのため表示しない．

### `visualize_archimedean_nets.m`（アルキメデス多面体版）

```mathematica
Get["/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/visualize_archimedean_nets.m"]
```

`archimedean_faceup_results.mx` を読み込み，
13種のアルキメデス多面体 × RS/RZ × with fallback の展開図を描画する．
`showR2 = True` にすると R2 (loxodrome) も表示される（デフォルト False）．
各多面体の見出しに RS/RZ/R2 の成功数・成功率サマリも表示する．

### `unfold4D.m`（4D展開図可視化，2026-05-08）

```mathematica
Get["/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/unfold4D.m"]
```

`ans4DGlobal.mx`（`peeling4Df4.m` の計算結果）と `unfold3DExport.m` を読み込み，
5種の正多胞体（5・8・16・24・120胞体）の 4D→3D 展開図を `Graphics3D` で描画する．

**固定 C1 の全 C2 バリアントを表示**：等変性により C1=3 を代表として固定（`fixedTop = 3`）し，
対応する全成功 C2 組み合わせを order でユニーク化して描画する．

**色付け**：青（C1）→緑（中間）→赤（末尾セル）のグラデーション（`apColor`）．

**5胞体の補正**：`ans4DGlobal.mx` に格納された `localVers` は f5.m 修正前の座標で計算されているため，
`recomputeLocalVers` で現在の f5.m から face-up 回転済み座標を再計算して使用する．
f5.m の修正後（2026-05-08）に `run4DPeelingGlobal.m` を再実行すれば，
`ans4DGlobal.mx` が更新されてこの補正は不要になる（ただしコードはそのままでも動作する）．

| 関数 | 役割 |
|------|------|
| `showNet4D[vers, faces, cells, order, sz]` | `unfoldTo3D` で展開し `Graphics3D` を返す |
| `apColor[k, nTotal]` | 青→緑→赤グラデーション |
| `getPairsForTop[nPoly, top]` | 固定 C1 の全成功ペアを取得，order で重複除去 |
| `showPolytopeNets[nPoly, poly, top, label, col, customVers]` | 1多胞体分を全 C2 バリアント表示 |
| `recomputeLocalVers[vers, faces, cells, top]` | 頂点座標更新後に face-up 座標を再計算 |
| `load4DData[n]` | `_4DData/f{n}.m` を読み込み Association を返す |

#### 4D 展開図におけるキラリティ正規化の検討（2026-05-12）

3D 多面体の 2D 展開図（`p3lUnfoldNet`）に対して実施したキラリティ正規化（外側から見て CW に統一）が
4D 多胞体の 3D 展開図にも必要かを検討した結果，**4D 展開図への適用は不要**と判断した．

理由：
- `p3lUnfoldNet` の問題は「2D 平面に射影するとき，外向き法線を +z に向けると表裏（手前/奥）が意味を持つ」という 2D 特有の問題
- `unfoldTo3D`（4D→3D 展開）は `Graphics3D` で表示するため視点自由であり，「表面が手前」という概念がない
- `unfoldProcrustes` が `det(R)=1` を強制するため，Procrustes 整列は常に正回転であり鏡映を混入しない
- SVD 射影基底の符号は任意だが，それは 3D 空間内での全体向きを決めるだけで，インタラクティブな回転で任意視点から観察できる Graphics3D では問題にならない

対応した修正（2026-05-12）：
- `unfold4D.m` の設定4行を `If[!ValueQ[...]]` で ValueQ 保護（`deduplicateGeo` の上書きバグ修正）

### `visualize_platonic_nets.m`（旧版・参照用）

R1/R2/R3 の3規則すべてを表示する旧版．アーカイブとして残してある．

### 主な関数

| 関数 | 役割 |
|------|------|
| `p3lUnfoldNet[vers, faces, order]` | ペーリング順序から2D展開図座標を生成 |
| `rotationMatrix3D[from, to]` | Rodrigues の公式による3D回転行列 |
| `peelingColor[k, nfTotal]` | 青（先頭）→緑（中間）→赤（末尾）のグラデーション色 |
| `netGraphicsOK[polys2D, nfTotal, sz]` | 成功展開図の描画（薄灰色背景） |
| `netGraphicsFail[polys2D, nPlaced, nfTotal, sz]` | 失敗展開図の描画（× マーク・配置面数/全面数ラベル） |
| `appendNetRows[...]` | 1規則分の with/no fallback 行を output に追記するヘルパー |

### `p3lUnfoldNet` のアルゴリズム

1. F1 の**外向き**法線が +z になるよう全頂点を回転（`Cross[e1,e2]` が内向きなら符号反転して外向きを保証）
2. F1 の xy 座標を2D座標として確定（面ローカル辞書 `fp[1, v]` に格納）
3. k=2..n：F(k-1) と Fk の共有辺を基準に，辺方向成分（t）と垂直距離（r）を3D座標から計算し，F(k-1) 重心と反対側に新頂点を配置（`fp[k, v]` に格納）
4. c1→c2 ベクトルが +x 方向になるよう最終回転

同じ頂点でも面ごとに独立した2D座標を持つ（`fp[k, v]` はモジュールローカルな DownValue）．

#### 表示規約

- **表示規約**：展開図を多面体の**外側から見たとき，螺旋が時計回り（CW）に見える**（右利きの人が左方向に剥く動作の視点）
- **根拠**：アルゴリズムが右半空間条件（det ≦ eps）を使うため，RS は 3D 空間で外側から見て CW の経路を直接選ぶ（F3 が y < 0 側）．**表示での y 反転は不要**
- **下半球 F1 との整合性**：`p3lAlignTopToZ` の `{1,−1,−1}` 変換（x 軸周り 180° 回転）が y を反転させるが，右半空間条件の下では F3 が引き続き y < 0 に収まるため上・下半球の F1 で同一レイアウトが得られる
- **実装ファイル**：`visualize_platonic_nets_v2.m` および `visualize_archimedean_nets.m` の `p3lUnfoldNet`（y 反転ブロック削除済み）

### 設定パラメータ（スクリプト冒頭で変更可）

| パラメータ | デフォルト | 説明 |
|-----------|:----------:|------|
| `netSize` | 120 | 各展開図の画像サイズ (px) |
| `maxShowOK` | `Infinity` | 成功展開図の表示上限 |
| `maxShowFail` | `Infinity` | 失敗展開図の表示上限 |
| `deduplicate` | `True` | `True`: 同じ order を1件に集約，`False`: 全ペア表示 |
| `deduplicateGeo` | `False` | `True`: 幾何学的重複除去（鏡像・回転等価を同一視）．`Get` 前に設定しないと上書きされないよう ValueQ 保護済み |
| `showStrict` | `True` | `True`: 厳密条件（eps < 0）の結果も表示（Platonic 版のみ） |
| `showR2` | `False` | `True`: R2 (loxodrome) も表示（Archimedean 版のみ） |
