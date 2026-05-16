# Opus 4.7 セッション作業ログ — 2026/05/14

`paper_draft.tex` の投稿準備に関する Opus モデルでの追加レビュー・修正作業の記録。

---

## 1. 全体評価（Opus による独立レビュー）

Sonnet モデルが行った修正後の論文を Opus 4.7 で再評価。**投稿水準にほぼ達しているが、用語と表記に複数の不整合が残存**と判定。

### Opus が指摘した主要問題

| # | 問題 | 重要度 |
|---|------|:------:|
| 1 | 「left-half-space」表記が 4 箇所残存（定義は right-half-space） | 重大 |
| 2 | Section 3.6 見出し「Spiral Path under RZ」が RS（Spiral rule）と用語衝突 | 重大 |
| 3 | 「Face-up」「cell-centroid-up」「3D-face-centroid-up」「face-centroid-up」が混在 | 中 |
| 4 | Algorithm 1・2 の `\textbf{exit}` / `\textbf{or}` 表記が if-else 構造として読みづらい | 中 |
| 5 | 引用キー `Yoshino2024arXiv` の年号が 2026 | 軽 |
| 6 | タイトルページ日付 2026/05/12 が古い | 軽 |
| 7 | CLAUDE.md の著者「Yoshino, Chaidee, Sawatdithep」が実際の論文（Yoshino, Chaidee）と不一致 | 軽 |

---

## 2. 適用した修正

### 2.1 CLAUDE.md と memory ファイルの著者修正

- `CLAUDE.md` line 5：「Yoshino, Chaidee, Sawatdithep」→「Yoshino, Chaidee」
- `memory/project_4d_peeling.md` line 7：同様の修正

### 2.2 用語の整合化（`paper_draft.tex`）

| 修正内容 | 詳細 |
|----------|------|
| left-half-space → right-half-space | 4 箇所（line 843, 1439, 1460, 1670） |
| Section 3.6 見出し変更 | `Spiral Path under RZ` → `Hexagonal-Band Path under RZ` |
| 文中の「spiral」表現修正（line 1010 周辺） | `a spiral around the peeling axis` → `a zonal band-by-band traversal around the peeling axis` |
| Archimedean Observations 内（line 798） | `produces a smooth spiral band` → `produces a smooth band-by-band traversal` |

### 2.3 4D 用語の正式定義と統一

**4D Setup 節（line 1029 付近）に正式定義を追加：**

```latex
The algorithm rotates the polytope so that the centroid of the starting
cell C_1 aligns with the +w direction; we call this orientation
**cell-centroid-up**, the 4D analogue of the 3D Face-up rotation.
...
An alternative orientation, **3D-face-centroid-up**, aligns +w with
the centroid of the shared 2-face between C_1 and a chosen neighbour;
this variant is examined in Section 4.3.
```

**他の箇所も統一：**
- `F_1` → `C_1`（4D Setup：cell に対する記号修正）
- 4D Remark：「Face-up orientation」→「Cell-centroid-up orientation」
- 「After Face-up alignment of C_1」→「After cell-centroid-up alignment of C_1」
- 600-cell 節：「face-centroid-up orientation」→「3D-face-centroid-up orientation」
- Algorithm 2 の `Rotate so that c_{C_1} aligns with +w` → `Rotate (cell-centroid-up) so that c_{C_1} aligns with +w`

### 2.4 Algorithm 1・2 の variant 表記改善

**修正前（混乱しやすい構造）：**

```
\Else
  \State \textbf{exit} (failure) \Comment{no-fallback variant}
  \State \textbf{or} (with-fallback variant):
  \If{$r = \text{RS}$}
  ...
```

**修正後（明確な if-elseif-else）：**

```
\If{R ≠ ∅}
  ... (rule-based selection)
\ElsIf{v = n}  \Comment{no-fallback variant}
  \State \Return (order, failure)
\Else  \Comment{with-fallback variant: R=∅, v=w}
  \If{r = RS}
  ...
```

**Algorithm caption も更新：**
- `Apple-Peel Unfolding (3D, rule r ∈ {RS, RZ})` → `Apple-Peel Unfolding (3D, rule r ∈ {RS, RZ}, variant v ∈ {w, n}). Variant w = with fallback; variant n = no fallback.`
- Algorithm 2 も同様

### 2.5 軽微な修正

- タイトル日付：`\date{2026/05/12}` → `\date{2026/05/14}`
- 引用キー：`Yoshino2024arXiv` → `Yoshino2026arXiv`（本文・bibliography 両方で統一）

---

## 3. 整合性チェック結果

すべての修正後、以下を機械的に検証：

- 未定義 `\ref` ：なし
- bibitem と cite の対応：全 10 件 1:1
- left-half-space・Sawatdithep・Yoshino2024arXiv の残存：なし

---

## 4. 論理的整合性の独立検証

「論理的矛盾がないか」について Opus が独立検証。**重大な論理矛盾は検出されず**。

### 検証で OK と確認した非自明な整合性

1. **8 胞体の Det=0 縮退**：centroid が ±e_i で det=0 が常に成立 → max-det が無意味化 → xy-cross-product フォールバックで 48/48 成功。Algorithm 2 line 1102–1106 と整合。
2. **Mirror symmetry の論理**：右螺旋（mirror）≡ 左螺旋（original）。amphichiral 11 種は equivariance、chiral 2 種は orbit decomposition で count 一致。
3. **C2 候補の w 値**：120 胞体の C1=3 で 12 個の C2 が band 8（w≈2.118）に位置 → k=2 で max-w タイ → xy-cross-product タイブレーク。
4. **9 帯構造**：1+12+20+12+30+12+20+12+1=120 算術整合。
5. **Strict ε<0 vs standard ε>0**：Remark 2.5 で明確に区別、両者とも equivariance 保持。

### 表現上の改善余地（論理矛盾ではない）

Conclusion 項目 2（"three solids"）と項目 3（"four solids"）が並んで紛らわしい可能性 → 修正実施。

---

## 5. Conclusion 部の明確化

**修正前：**
```
2. RZ achieves 100% on three Archimedean solids
   (Truncated Octahedron, Truncated Icosahedron, Truncated Cuboctahedron);
   seven of thirteen remain Impossible under both rules.
3. RZ outperforms RS on four solids: ...
```

**修正後：**
```
2. RZ achieves 100% on three of thirteen Archimedean
   solids (Truncated Octahedron, Truncated Icosahedron,
   Truncated Cuboctahedron), and 66.7% on the Truncated Icosidodecahedron;
   seven of thirteen remain Impossible under both rules.
3. RZ outperforms RS on these four solids: ...
```

「three」と「four」の関係が「3 種の 100% + 1 種の 66.7% = these four」として論理的に接続。

---

## 6. 投稿先ジャーナル推奨（Opus 提案）

### 第一候補：**Computational Geometry: Theory and Applications（CGTA, Elsevier）**

- 引用文献 Devadoss & Harvey (2022) が掲載された誌
- 内容適合度が最高、リジェクトリスクが相対的に低い
- Scopus 収録確実、評価サイクルが業界標準

### 第二候補：**Discrete & Computational Geometry（DCG, Springer）**

- 分野最高峰、引用 O'Rourke (2015) や Aronov & O'Rourke (1992) が掲載
- 等変性証明（Remark 3.1）が評価される可能性
- ただし査読が非常に厳しい

### 第三候補（隠れた良候補）：**Discrete Mathematics（Elsevier）**

- 引用 Buekenhout & Parker (1998) が掲載された誌
- 4 次元正多胞体の net 列挙という本論文の中核テーマと直接系統
- 組み合わせ論的読者が 4D 結果を歓迎する可能性

### オープンアクセス候補：**Symmetry（MDPI）**

- 対称性・等変性の議論が中心の本論文に適合
- 査読が早い（1〜2 ヶ月）
- 掲載料 CHF 2,000〜

### 回避推奨

- 一般組合せ論誌（European Journal of Combinatorics 等）：アルゴリズム面が薄まる
- Symmetry: Art and Science（Kaino 2019 の掲載誌）：Scopus 非収録の可能性が高い
- Mathematics 等の一般数学誌：査読の質が読めない

### 最終推奨

**第一投稿先：CGTA**

理由：
1. 内容適合度が最高、リジェクトリスクが相対的に低い
2. 既往研究との系統（Devadoss2022）があり査読者が背景を理解しやすい
3. Scopus 収録確実

**DCG に挑戦する場合**：120 胞体の自己交差解析と等変性証明をさらに深め、Discussion を厚くしてから投稿することを推奨。

---

## 7. 投稿準備状況

すべての修正後、論文は以下の状態：

- **投稿水準に達している**
- 用語と表記が一貫
- 数値・表・本文の主張が高度に一貫
- 等変性論証・mirror symmetry 分析・Platonic/Archimedean 分類・4D 結果がすべてデータによって支持

次のステップは投稿先の最終決定とジャーナル指定のフォーマット調整。
