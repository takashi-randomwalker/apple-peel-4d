# 4次元 Apple-Peel 展開：現状まとめと再開手順

**最終更新：** 2026-04-03
**ステータス：** 一時中断中（3D規則比較の検討を優先）

---

## 1. プロジェクト概要

3次元多面体の Apple-Peel 展開（Yoshino, Chaidee, Sawatdithep）を
**4次元正多胞体**に拡張するプロジェクト．

- 4次元版の「左側条件」と選択規則を定式化し，Mathematica で実装
- 全6種の正多胞体を対象に展開可能性（Perfect / Possible / Impossible）を分類
- 展開図を SVD 射影 + Procrustes 整列で3次元可視化し，STL 出力

---

## 2. 実装ファイル一覧

| ファイル | 役割 | 状態 |
|----------|------|------|
| `peeling4Df.m` | 4Dペーリング本体（現行版，xy 2成分左条件 + 規則3） | 完成・使用中 |
| `peeling4Df3.m` | 改良版（xyz 行列式左条件 + 規則3）| **新規作成済み・未テスト** |
| `run4DPeeling.m` | 全6胞体の実行スクリプト（`peeling4Df.m` 使用） | 完成 |
| `run4DPeelingR3.m` | 改良版の実行スクリプト（`peeling4Df3.m` 使用） | **新規作成済み・未実行** |
| `peeling1back.m` | 1段バックトラック拡張 | 完成・未比較 |
| `unfold3DExport.m` | SVD射影・Procrustes整列・SAT判定・STL出力 | 完成 |
| `_4DData/f5.m` 〜 `f600.m` | 各正多胞体のデータ（頂点・面・セル） | 完成 |
| `ans4D.mx` | 現行版の計算結果（Mathematica バイナリ） | 保存済み |
| `STLs/` | 570枚のSTLファイル | 保存済み |
| `summary.tex` / `.pdf` | アルゴリズム詳細ドキュメント | 最新版（下記参照） |

---

## 3. 計算結果（`peeling4Df.m` 現行版）

左条件：$c_{k,y} c'_x - c_{k,x} c'_y \ge 0$（xy 2成分版）
選択規則：規則3（左候補から max $w$，フォールバック min $w$）

| 多胞体 | セル数 | 隣接数/セル | $(C_1,C_2)$ 総数 | ユニーク成功数 | 分類 |
|--------|:------:|:-----------:|:----------------:|:-------------:|------|
| 5胞体  | 5  | 4 | 20    | 20  | **Perfect** |
| 8胞体  | 8  | 6 | 48    | 48  | **Perfect** |
| 16胞体 | 16 | 4 | 64    | 26  | Possible |
| 24胞体 | 24 | 8 | 192   | 119 | Possible |
| 120胞体 | 120 | 12 | 1,440 | 357 | Possible |
| 600胞体 | 600 | 4 | 2,400 | — | **未完了** |

---

## 4. `summary.tex` の構成（最新版）

| 節 | タイトル | 内容 |
|----|----------|------|
| §1 | 背景と目的 | Apple-Peel展開の説明，対象正多胞体 |
| §2 | 3次元アルゴリズムの概要 | アルゴリズム1（擬似コード付き） |
| §3 | 4次元への拡張 | w軸整列，隣接判定，4D回転 |
| **§4** | **左側条件の再検討と4次元への厳密な拡張** | **今回新規追加（下記参照）** |
| §5 | `peeling4Df` の実装 | アルゴリズム2（擬似コード付き） |
| §6 | ファイル構成と使い方 | |
| §7 | 計算結果 | 分類表，成功率 |
| §8 | 展開図の3次元可視化と検証 | SVD・Procrustes・SAT・STL数 |
| §9 | 1段バックトラック拡張 | アルゴリズム3（擬似コード付き） |
| §10 | 今後の課題 | 下記と同内容 |
| 付録 | `peeling4Df` 完全コード | |

### §4 新規追加内容（重要）

`summary.tex` の §4「左側条件の再検討と4次元への厳密な拡張」に以下を追加済み：

- **§4.1** Darboux フレームを用いた3D左条件の正確な再定式化
  　→ $\det(\mathbf{c}_{k-1}, \mathbf{c}_k, \mathbf{c}_j) \ge 0$ が $\hat{L}_k \cdot \mathbf{c}_j \ge 0$ と完全等価であることを証明
- **§4.2** $\mathbb{R}^4$ における困難：直交補空間が2次元になるため「左」が一意に定まらない
- **§4.3** 定義4.1：直前ステップ参照による4次元左条件
  　→ $\det(\tilde{\mathbf{c}}_{k-1}, \tilde{\mathbf{c}}_k, \tilde{\mathbf{c}}_j) \ge 0$（$xyz$ 3成分）
  　→ $k=1$ では直前セルが未定義のため xy 2成分条件を使用（初期条件の不定性）
- **§4.4** 3つの選択規則の定義：規則1（max $\varphi$），規則2（min $\varphi$ / loxodrome），規則3（max $w$）
- **§4.5** loxodrome vs 緯度優先の幾何学的比較
- **§4.6** Kaino (2019) との対応（規則3 に相当）
- **§4.7** 2段階選択の提案：規則3 優先 → 規則2 補完（表付き）

---

## 5. 再開時の最優先タスク

### 5.1 `peeling4Df3.m` のテスト（最重要）

新しい行列式左条件を実装した `peeling4Df3.m` / `run4DPeelingR3.m` が未実行．
まず小規模な5胞体・8胞体で旧版と結果を比較する：

```mathematica
SetDirectory["...260324Peeling4D"];
Get["peeling4Df3.m"];
(* 5胞体で確認 *)
raw5 = Get["_4DData/f5.m"];
{vers, _, _, faces, _, cells} = raw5;
res_old = peeling4Df[N[vers],  faces, cells, 1];
res_new = peeling4Df3[N[vers], faces, cells, 1];
(* 結果比較 *)
{Last /@ res_old, Last /@ res_new}
```

### 5.2 全6胞体の再計算

```mathematica
Get["peeling4Df3.m"];
Get["run4DPeelingR3.m"];
(* results4DR3 に結果が格納される *)
```

### 5.3 その他の課題（優先度順）

| 優先度 | タスク |
|:------:|--------|
| 高 | `peeling4Df3.m` のテスト（旧版との比較） |
| 高 | 600胞体の計算完了（計算量大，分割実行） |
| 中 | STL 570枚の重なりなし（valid net）比率集計 |
| 中 | `peeling1back.m` との成功率比較 |
| 中 | 2段階選択（規則3 → 規則2）の実装と評価 |
| 低 | Kaino (2019) との比較 |
| 低 | 左利きペーリングとの対称性検討 |

---

## 6. 4次元左条件の数学的背景（要点）

$$
\text{3Dの左条件（面重心 } \in \mathbb{R}^3 \text{）：}
\quad \det(\mathbf{c}_{k-1}, \mathbf{c}_k, \mathbf{c}_j) \ge 0
$$

$$
\text{4Dの左条件（xyz 射影 } \tilde{\mathbf{c}} \in \mathbb{R}^3 \text{）：}
\quad \det(\tilde{\mathbf{c}}_{k-1}, \tilde{\mathbf{c}}_k, \tilde{\mathbf{c}}_j) \ge 0
$$

- 現行実装（`peeling4Df.m`）は $c_{k,y}c'_x - c_{k,x}c'_y \ge 0$（xy 2成分のみ，$k$ によらず同じ式）
- 改良版（`peeling4Df3.m`）は $k \ge 2$ で上記3成分行列式，$k=1$ は xy 2成分（直前セル未定義のため）

---

## 7. 関連ファイル・参考情報

- **`results_3D_rule_comparison.md`** — 3D正多面体での規則2 vs 規則3 比較実験結果
  （正12面体で規則2が常に6面目で詰まる「6面の壁」を確認）
- **`peeling3DLoxo.m`** — 3D版規則2・規則3の Mathematica 実装
- **`summary.pdf`** — 最新の全アルゴリズム記述（PDF）
