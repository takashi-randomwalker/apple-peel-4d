# Opus 4.7 レビュー（2026-05-16）— `paper_draft.tex`

## 結論

- **arXiv 投稿は即可**
- **CGTA（第1候補）は submittable，ただし major revision を覚悟すべき水準**
- **DCG レベルを狙うなら追加の形式化（命題化・600胞体の structural argument）が必要**

---

## 投稿水準に達している点

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

---

## 査読で指摘される可能性が高い点（major revision 級）

### 1. 等変性結果が「Remark」止まり

Section 2 の等変性議論は数式は揃っているが正式な命題化されていない。CGTA・DCG では：

```latex
\begin{proposition}[Equivariance]
Let G ⊂ SO(3) be the symmetry group of a face-transitive polyhedron.
Algorithm 1 with rule r ∈ {RS, RZ} is equivariant under G, so the success
indicator depends only on the G-orbit of (F_1, F_2).
\end{proposition}
\begin{proof} ... \end{proof}
```

の形にすべき。これがないと「informal で computational」と評される。

### 2. 600-cell の "icosahedral bottleneck" が経験則止まり

5 通りの停止ステップが固定（146, 150, 276, 279, 284）という観察はあるが，それが正二十面体対称に由来することの構造的説明が薄い。少なくとも：

- 各停止ステップで Det 条件が「全候補を排除する」具体例を1つ示す
- セル隣接グラフの構造的特徴（接続数4 + I_h 対称）を補題化する

のいずれかは必要。

### 3. コード/データ可用性ステートメントの欠如

Mathematica 実装（`peeling3DLoxo.m`, `peeling4Df4.m`, `_4DData/`）の公開予定を明記すべき。現状記載なし。近年の computational geometry 系誌は事実上必須：

```latex
\section*{Code and Data Availability}
The Mathematica implementation and computed result archives are available at
https://github.com/.../apple-peel-4d (or Zenodo DOI).
```

### 4. 参考文献が9件と寡少

最低でも以下は追加検討：

- Schlickenrieder, *Nets of Polyhedra* (TU Berlin thesis, 1997) — spiral unfolding の起源
- Pak, *Lectures on Discrete and Polyhedral Geometry* (2010)
- Bern, Demaine, Eppstein, et al. *Ununfoldable polyhedra* (1999)
- 4D unfolding 関連の Coxeter, Towle など歴史的言及

### 5. 「なぜこの2規則か」の動機付けが弱い

Section 2.2 で RS と RZ が天下りに定義されている。「他の plausible rule（例：min angle, nearest neighbor）と比べてなぜこの2つを選んだか」を1段落追加すべき。

---

## minor だが対応すれば質が上がる

- **Section 5（Computational Examples）が Section 3, 4 と一部重複** — 個別深掘りに整理した方がよい。特に Section 5.1（Cross-Dimensional Summary）の表は Section 3, 4 の結果の集約のみ。
- **Conclusion の Future Work が短い** — 「two-stage selection の計算結果を将来公開」「random polyhedra での評価予定」など，より具体的に。
- **Abstract に "we classify each solid as Perfect/Possible/Impossible" とあるが，この分類が新しい貢献の一部であることを明記**するとよい。
- **Acknowledgements に共著者の funding 情報があるが ORCID, conflict-of-interest 宣言**は誌指定で必要になる場合あり。

---

## 推奨アクション

| 投稿先 | 必要な追加作業 | 期間目安 |
|--------|---------------|---------|
| **arXiv（即時）** | なし | 即日 |
| **CGTA**（第1候補） | 上記 1–3 を対応（Proposition 化 + 600-cell の補強 + Code availability） | 1–2 週間 |
| **DCG**（高目標） | 上記すべて + 600-cell の Impossibility に近い structural argument | 1–2 ヶ月 |
| **J. Comput. Geom.** | 1, 3 のみ対応で十分 | 数日 |

**推奨フロー**：arXiv に先に置いて DOI を確保し，並行して CGTA 向けに 1, 3 だけ最低限対応してから submit する流れが効率的。3 (code availability) は GitHub レポジトリを準備するだけなので 1 日で済む。1 (Proposition 化) は既存の議論を整形するだけなので半日〜1日。これで査読初稿提出は可能水準になる。

---

## このセッション（2026-05-16）で行った修正

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

---

## レビュー方法

- 主モデル：Claude Opus 4.7（`/model` で切替）
- 補助モデル：Claude Sonnet 4.6（事実確認のサブエージェント）
- 対象ファイル：`/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/paper_draft.tex`（1,710 行，30 ページ PDF）
- レビュー観点：数値正確性，用語一貫性，図表参照，文法・表現，構造的問題，R2 残存チェック，投稿水準判定
