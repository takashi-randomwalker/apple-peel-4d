# HISTORY.md — 4D Apple-Peel 展開プロジェクト 履歴・廃止実装・詳細ログ

このファイルは `CLAUDE.md` から退避した詳細記録です。
現行版の仕様・最新結果は `CLAUDE.md` を参照してください。

---

## 廃止された実装

### `peeling4Df.m`（xy 2成分外積，旧版）

左条件：`c_{k,y} * c'_x - c_{k,x} * c'_y >= 0`（xy 2成分外積）

| 多胞体 | セル数 | 隣接数/セル | (C1,C2)総数 | ユニーク成功数 | 分類 |
|--------|:------:|:-----------:|:-----------:|:-------------:|------|
| 5胞体  | 5      | 4           | 20          | 20            | **Perfect** |
| 8胞体  | 8      | 6           | 48          | 48            | **Perfect** |
| 16胞体 | 16     | 4           | 64          | 26            | Possible |
| 24胞体 | 24     | 8           | 192         | 119           | Possible |
| 120胞体 | 120   | 12          | 1,440       | 357           | Possible |
| 600胞体 | 600   | 4           | 2,400       | —             | **未完了** |

### `peeling4Df3.m`（xyz ローカル k-1 参照）

左条件：`Det[{c_{k-1}, c_k, c_j}]_{xyz} >= 0`（k=1 は xy 2成分にフォールバック）
作成済みだが未テストのままアーカイブ。

### `peeling4Df4.m` v3（タイブレーク前版，2026-05-10）

計算スクリプト：`run4DGlobal_xy2.m`，結果：`ans4DGlobal_v3.mx`

| 多胞体 | (C1,C2)総数 | True | Unique | 分類 |
|--------|:-----------:|-----:|-------:|------|
| 5胞体  | 20          | 20   | 20     | Perfect |
| 8胞体  | 48          | 48   | 48     | Perfect |
| 16胞体 | 64          | 25   | 25     | Possible |
| 24胞体 | 192         | 192  | 192    | Perfect |
| 120胞体 | 1,440      | 1,250| 1,250  | Possible |
| 600胞体 | 2,400      | 0    | 0      | Impossible |

v3 → v4 の変化：16胞体 25→20，8胞体 geo-unique 5→1，**120胞体 1250→1440（Perfect）**。
これは k≥3 の w タイブレーク時に Det が全候補で ≈0 の場合に xy クロス積へフォールバックするよう修正したため。

### `peeling1back.m`（旧版1段バックトラック）

xy ローカル参照で実装。`peeling4Df4_1back.m`（Det グローバル参照版）に置き換えられ使用非推奨。

### `STLs/`（570枚），`STLs_v3/`（285枚）

- `STLs/`：旧版，重なり未チェック
- `STLs_v3/`：v3 結果からの valid net（重なりなし保証）
- 現行は `STLs_v4/`（280枚，2026-05-12）

### `ans4DGlobal_v2.mx`

タイブレーク後の計算結果。DumpSave 失敗で 0 バイト。
DumpSave の既知の問題：`results4D_RZ[n] = <|...|>` のように Association の DownValue を DumpSave すると 0 バイトになる。`results4DGlobal` という DownValue 形式は問題なく動作。

---

## 1段バックトラック結果（`peeling4Df4_1back.m` + RZ，2026-05-11）

計算スクリプト：`run4DBacktrack.m`，結果：`ans4DGlobal1back.mx`

| 多胞体 | (C1,C2)総数 | 1back True | greedy True | 差 |
|--------|:-----------:|:----------:|:-----------:|:---:|
| 5胞体  | 20          | 20         | 20          | 0   |
| 8胞体  | 48          | 48         | 48          | 0   |
| 16胞体 | 64          | 25         | 25          | 0   |
| 24胞体 | 192         | 192        | 192         | 0   |
| 120胞体 | 1,440      | 1,250      | 1,250       | 0   |
| 600胞体 | 2,400      | 0          | 0           | 0   |

- **1段バックトラックはグリーディ版と完全に同一**
- 600胞体は Impossible のまま（構造的ボトルネックはローカル BT で解消不可）
- アルゴリズム：Det グローバル参照 + RZ；justBT フラグで連続バックトラック禁止
- 計算時間：600胞体 516 秒（8.6 分）

> 注：上の数値は v3 比較時点のもの。v4 fix 後の現行値は 120胞体 1,440 で，1back も同様に 1,440 となる。

---

## 600胞体の終了ステップ分布（RZ，2026-05-10）

計算スクリプト：`analyze_600cell_stuck.m`（全2,400ペアを直接計算）

| 終了ステップ | ペア数 | 割合 |
|:-----------:|-------:|-----:|
| 146         | 742    | 30.9% |
| 150         | 825    | 34.4% |
| 276         | 404    | 16.8% |
| 279         | 323    | 13.5% |
| 284         | 106    |  4.4% |

全 2,400 ペアがわずか 5 通りの終了ステップで停止 — 正二十面体対称に由来する構造的ボトルネック。

## 600胞体の cell-centroid-up w 分布（2026-05-12）

計算スクリプト：`analyze_600cell_wdist.m`，C1=1 で face-up 後の全 600 セル重心 w 座標を解析

### cell-centroid-up の層構造（31 バンド）

| 層（上から） | バンド数累積 | w 値 | 備考 |
|:------------:|:-----------:|:-----:|------|
| band 9 終端  | 137 | +0.926 | |
| band 10 終端 | 161 | +0.749 | |
| band 15 終端 | 273 | +0.177 | |
| band 16（赤道帯） | 327 | ≈0 | **54 セル（最大バンド）** |

C2 候補（4 個）の w 座標は完全一致（spread < 10⁻¹⁵）— 600胞体でも C2 w 等値定理が成立。

### 停止ステップと層境界の対応

| 停止ステップ | 最近接（cell-up） | 差 | 最近接（Kaino） | 差 |
|:-----------:|:----------------:|:--:|:---------------:|:--:|
| 146 | 137（band 9 終端） | 9 | 130（layer 4） | 16 |
| 150 | 161（band 10 終端） | 11 | 130（layer 4） | 20 |
| 276 | 273（band 15 終端） | **3** | 270（layer 6） | 6 |
| 279 | 273（band 15 終端） | 6 | 270（layer 6） | 9 |
| 284 | 273（band 15 終端） | 11 | 270（layer 6） | 14 |

- Kaino の vertex-first 層境界との対応という当初仮説は**不支持**（差が 2 倍程度大きい）
- cell-centroid-up 自身の層境界との対応の方が有意に良い
- **ステップ 276-284 ≈ 累積 273**：w≈0.177 の band 15 末尾 → 直後に赤道帯（w≈0，54 セル）への遷移で詰まる構造的バリア
- **ステップ 146-150 ≈ 累積 137-161**：w≈0.926〜0.749 の中間帯での別の構造的バリア
- 接続数の制約：600 胞体は隣接数 4（120 胞体は 12）→ Det 条件がさらに絞ると行き詰まりやすい

---

## 3D 計算結果の詳細

### 3D 正多面体（Face-up あり，fallback あり/なし）

| 多面体 | ペア数 | R1(w) | R1(n) | R2(w) | R2(n) | R3(w) | R3(n) |
|--------|:------:|------:|------:|------:|------:|------:|------:|
| Tetrahedron | 12 | 100% | 100% | 100% | 100% | 100% | 100% |
| Cube | 24 | 100% | 100% | 100% | 0% | 100% | 100% |
| Octahedron | 24 | 100% | 100% | 100% | 0% | 100% | 100% |
| Dodecahedron | 60 | 100% | 100% | 100% | 0% | 100% | 100% |
| Icosahedron | 60 | 100% | 100% | 0% | 0% | 100% | 100% |

w = with fallback，n = no fallback。アルゴリズム：大域的 c_1 参照・eps 閾値・単一候補優先なし（2026-05-06 版）。
全結果が 0% か 100% — 面推移性による等変性理論と一致。

### 3D 正多面体（厳密条件 eps < 0, thresh = +10⁻¹⁰）

計算スクリプト：`run_platonic_strict.m`，結果：`dataPlatonic_strict.mx`（2026-05-07）

| 多面体 | ペア数 | RS(w) | RS(n) | RZ(w) | RZ(n) |
|--------|:------:|------:|------:|------:|------:|
| Tetrahedron | 12 | 100% | 100% | 100% | 100% |
| Cube | 24 | 100% | 0% | 100% | 0% |
| Octahedron | 24 | 100% | 0% | 100% | 0% |
| Dodecahedron | 60 | 100% | 0% | 100% | 0% |
| Icosahedron | 60 | 100% | 0% | 100% | 0% |

- with fallback は常に 100%（「真南面」が fallback 経由で選ばれる）
- no fallback は Tetrahedron のみ 100%（正四面体には「真南面」が存在しない）
- RS と RZ が完全一致（正多面体の緯度帯構造による）

### 3D アルキメデス（詳細，2026-05-12 Det ベース）

計算スクリプト：`run_archimedean_faceup.m`（fallback あり），`run_archimedean_nofallback.m`（fallback なし）

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

- **7 種（54%）が RS・RZ とも 0%**
- **RZ > RS の 4 種**：TruncatedOctahedron, TruncatedCuboctahedron, TruncatedIcosahedron, TruncatedIcosidodecahedron
- **RS > RZ の 1 種**：TruncatedTetrahedron
- **RS = RZ は 8/13 種（いずれも 0%）**
- RZ 結果は Det ベース選択導入前後で完全不変（min-φ ≡ min-Det）
- RS 結果は大幅変化（旧版で 100% だった TruncatedCuboctahedron/TruncatedIcosahedron/TruncatedIcosidodecahedron が Det ベースで 0%）
- フォールバック効果：SnubCube のみ（RS: 24→0，RZ: 48→24）

---

## Face-up 方式の比較実験（120胞体 C1=3）

スクリプト：`analyze_120cell_c2types.m`，`test_faceup_pair.m`，`test_faceup_allC2.m`，`gen_120cell_faceup_nets.m`

### Face-up 方式の比較

| 方式 | C2 w-bands | spread |
|------|:---------:|-------:|
| cell-centroid-up（現行）| 1（全同値）| ~10⁻¹⁵ |
| vertex-up | 4（3+3+3+3）| 0.926 |
| 3D-face-centroid-up | 4（1+5+5+1）| 0.951 |

**3D-face-centroid-up**：C1 と C2 の共有ペンタゴン面の重心を +w に揃える。
3D の face-up（開始面の重心を +z）の直接類推であり，正十二面体の 1+5+5+1 緯度帯構造を復元。

### 全 12 C2 での比較（v4 fix 適用後，2026-05-13）

| C2 | A: cell-up | A-ov | B: face-up | B-ov | B−A |
|----|:----------:|:----:|:----------:|:----:|:---:|
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
| **合計** | 12/12 | **308** | 7/12 | **130** | — |

- B 成功（7 ケース）：C2 = 1, 2, 18, 19, 51, 97, 100
- B 失敗（5 ケース）：C2 = 4, 10, 15, 17, 115
- 旧テスト（2026-05-11，10 通り）は `peel1Pair` に v4 fix 未適用で C2=19/51 が漏れ，C2=115 は OK だった
- `STLs_120cell_faceup/`：成功 7 ケース × cell-up/face-up = 14 STL

### 600 胞体への 3D-face-centroid-up 適用（2026-05-11）

スクリプト：`test_600cell_faceup.m`（全 2,400 ペア，227 秒）

| 方式 | 成功数 | 停止ステップ種類 | 主な停止 |
|------|:------:|:--------------:|---------|
| cell-centroid-up（現行）| 0/2400 | 5 通り | 146, 150, 276, 279, 284 |
| 3D-face-centroid-up     | 0/2400 | 25 通り | 88〜201 に分散 |

- 停止ステップの多様化（5 → 25）はボトルネックパターンの変化を示すが Impossible 不変
- 最大停止が 284 → 201 に短縮（より早く詰まる）
- 主な停止点：step 99（320 件, 13.3%），step 103（319 件, 13.3%）

### 600 胞体への vertex-up 適用（2026-05-12）

スクリプト：`test_600cell_vertexup.m`（C1=1，4 頂点 × 4 C2 = 16 ケース）

| 方式 | 成功数 | 停止ステップ範囲 |
|------|:------:|:----------------:|
| cell-centroid-up | 0/4 | 146, 284（2 パターン） |
| vertex-up | 0/16 | 80〜206（14 種類） |

- vertex-up は方向の多様性をもたらすが w-based 道案内が不正確になり逆効果
- 多くのケースで cell-centroid-up より**早く**詰まる

### Up 方向変更の網羅的まとめ

| up 方向 | 規模 | 成功数 | 停止パターン |
|---------|:----:|:------:|------------|
| cell-centroid-up | 2,400 ペア | 0 | 5 固定（正二十面体対称） |
| 3D-face-centroid-up | 2,400 ペア | 0 | 25 種類に分散 |
| vertex-up | 16 ケース | 0 | 14 種類に分散，平均早く詰まる |

**結論**：up 方向の変更は投影構造を変えるが，接続数 4 という貧しいグラフ構造の中で Det 条件が全候補を封鎖する根本障壁を解消できない。

---

## 120 胞体の自己交差構造解析（2026-05-14）

スクリプト：`analyze_120cell_overlap_bands.m`，`analyze_120cell_band_backtrack.m`

### 帯間 vs 帯内の重なり分布

C1=3 の全 12 C2 × 各 308 重なりペアを解析：
- **同一帯内：88 ペア（28.6%）**
- **異なる帯間：220 ペア（71.4%）**

異帯重なりはほぼすべて隣接帯間（帯 k ↔ 帯 k+1）。自己交差は帯境界で各セルが隣接帯に「はみ出す」ことで発生。

### 帯境界代替エントリー検証

仮説：帯遷移ステップで RZ が選ぶセル以外を選べば重なりを回避できるか？

C1=3，C2=1 で実行：
- **全帯遷移ステップで「代替なし」（隣接する次帯セルが 1 個のみ）**
- RZ は帯境界を越える際に選択の余地がない
- RZ の帯走査は単純な「帯を順番に一周」ではなく，帯 4↔5，5↔6，6↔7，7↔8 の間で何度もジグザグ
- 帯遷移は約 60 回（全 119 ステップ中）で各遷移に唯一の通路

**結論**：帯境界での局所的セル置換で重なりは回避できない。解消にはグローバルな最適化が必要（グリーディ系では困難）。

---

## 数値誤差問題の解決（2026-05-06）

旧実装では Dodecahedron R3(with fallback) = 53.3% など，等変性理論（0% or 100%）に矛盾する結果が生じていた。

### 原因
1. 左条件の参照点が局所的（c_{k-1}）だったため，螺旋方向が各ステップで変化し対称性が崩れていた
2. 正十二面体の黄金比座標を `N[...]` で有限精度化すると，理論上 0 となるはずの行列式値が ~10⁻¹⁶ の誤差を持ち，符号判定がペアによって逆転していた
3. 単一候補優先が左条件をバイパスし，大域螺旋保証を破っていた

### 対処
- 左条件を `Det[{cc[[top]], cc[[last]], cc[[j]]}] >= -eps`（c_1 大域固定参照・`eps = 10⁻¹⁰`）に変更
- 単一候補優先を廃止（全候補数に関わらず左条件を適用）

### 結果
全 5 種の正多面体で全結果が 0% か 100% になり，等変性理論と整合。

---

## 対称性に関する理論的考察

### 正多面体における面推移性とアルゴリズムの等変性

正多面体は面推移的（face-transitive）であり，任意の 2 面を写す回転対称が存在する。

**結論：アルゴリズムは対称群 G の下で等変であり，成功率は 0% か 100% のいずれかでなければならない。**

証明の概略：σ ∈ G が (F1, F2) → (F1', F2') を与えるとき，

1. Face-up 回転後の重心：`cc_{F1'}[σ(j)] = A · cc_{F1}[j]`（A は proper rotation，det = 1）
2. 左半空間条件：`Det[{A·a, A·b, A·c}] = det(A)·Det[{a,b,c}] = Det[{a,b,c}]` → 左候補集合は σ で置換されるだけ
3. R1/R2（Det 選択）：`Det[{A·c_1, A·c_k, A·c_j}] = Det[{c_1,c_k,c_j}]` → min/max Det の順序は σ で不変
4. R3（max z）：z 座標は σ で面ラベルが置換されるだけ

したがって (F1, F2) で成功 ⟺ (F1', F2') で成功。面推移的正多面体では全ペアが等価なので成功率は 0% か 100% のみ。

---

## 4D 展開図におけるキラリティ正規化の検討（2026-05-12）

3D 多面体の 2D 展開図（`p3lUnfoldNet`）に対して実施したキラリティ正規化（外側から見て CW に統一）が 4D 多胞体の 3D 展開図にも必要かを検討した結果，**4D 展開図への適用は不要**と判断。

理由：
- `p3lUnfoldNet` の問題は「2D 平面に射影するとき外向き法線を +z に向けると表裏（手前/奥）が意味を持つ」という 2D 特有の問題
- `unfoldTo3D`（4D→3D 展開）は `Graphics3D` で表示するため視点自由
- `unfoldProcrustes` が `det(R)=1` を強制するため，Procrustes 整列は常に正回転で鏡映を混入しない
- SVD 射影基底の符号は任意だが，それは 3D 空間内での全体向きを決めるだけ

対応した修正：`unfold4D.m` の設定 4 行を `If[!ValueQ[...]]` で ValueQ 保護（`deduplicateGeo` の上書きバグ修正）。

---

## 120 胞体 RZ nets の幾何学的多様性（C1=3，2026-05-14 詳細）

スクリプト：`analyze_c2_xyz.m`

### Geo-unique net 数と軌道構造

| 軌道 | C2 メンバー | サイズ |
|------|:----------:|:------:|
| orbit 1 | 1, 15, 17, 97 | 4 |
| orbit 2 | 4, 18, 51 | 3 |
| singleton × 5 | 2, 10, 19, 100, 115 | 各 1 |

### C2 重心の正二十面体配置

Cell-centroid-up 後の 12 個の C2 セル重心（xyz 部分）は，黄金比 φ ≈ 1.618 を使った循環座標
$(0, \pm\varphi/2, \pm\varphi^2/2)$ の全置換により正十二面体の面配置に対応した正二十面体を形成（辺長 φ，各 r_xyz ≈ 1.539）。

### 対蹠ペア

12 個の C2 は 6 対の対蹠ペア：{1,18}, {2,17}, {4,19}, {10,115}, {15,100}, {51,97}。
対蹠ペアは k=2 タイブレーカー（xy クロス積スコア）の符号が逆転するため，同じ軌道に属さない。

### 普遍的内部コア（k=1..13）

7 種の net すべて，最初の 13 セル（C1 + 12 個の直接隣接セル）について r_3D が完全一致：

```
r_3D[k=1..13] = {0.00, 1.70, 1.79, 1.80, 1.81, 1.82, 1.81, 1.96, 1.99, 2.04, 2.06, 2.02, 2.09}
```

k=14 以降で分岐。

### 螺旋構造パターン

| C2 | 軌道 | 旋回数 | 方向 | z 範囲 | パターン |
|----|:----:|:------:|:----:|:------:|:--------:|
|  1 | {1,15,17,97} | +1.07 | CCW | −3〜+26 | A（CCW 上昇） |
| 10 | singleton | +0.85 | CCW | −5〜+8  | A（CCW 上昇） |
| 19 | singleton | +1.08 | CCW | −4〜+20 | A（CCW 上昇） |
|100 | singleton | +0.75 | CCW | −3〜+29 | A（CCW 上昇） |
|115 | singleton | +0.82 | CCW | −13〜+6 | A（CCW 下降） |
|  2 | singleton | −1.14 | CW  | −7〜+4  | B（CW 反転） |
|  4 | {4,18,51} | −0.07 | ≈0  | −4〜+23 | C（柱状） |

- 型 A（5/7 net）：RZ の max-xy-cross-product が CCW を選好することと整合
- 型 B（1/7，C2=2）：CW 反転型。C2=17 と対蹠
- 型 C（1/7，C2=4 代表）：旋回ほぼゼロ，z 方向に伸びる柱状

---

## 完了タスク履歴（時系列，参考）

| 日付 | タスク |
|:----:|--------|
| 2026-05-06 | `peeling3DLoxo.m` 大域 c_1 参照・eps・単一候補優先廃止 |
| 2026-05-07 | 厳密条件計算；アルキメデス再計算；SnubCube バグ修正 |
| 2026-05-08 | `peeling4Df4.m` 新規作成（c1–c2 平面グローバル参照）；`run4DPeelingGlobal.m` で全6胞体計算 |
| 2026-05-10 | RS を `peeling4Df4.m` に実装（`run4DGlobal_RSRZ.m`）：4D では RZ が RS を大きく上回る；600胞体の終了ステップ頻度分布解析；k=2 タイブレーク幾何学的基準（xy クロス積）；`run4DGlobal_xy2.m` で全 C1 計算 → `ans4DGlobal_v3.mx`；アルキメデスフォールバックなし計算 |
| 2026-05-11 | 1段バックトラック（`peeling4Df4_1back.m`，グリーディ版と完全に同一）；3D-face-centroid-up を 600胞体に適用（Impossible 不変）；STL valid net 比率集計（`STLs_v3/` 285枚） |
| 2026-05-12 | k≥3 Det=0 縮退で xy クロス積フォールバックを実装 → 120胞体 1250→1440（Perfect）；4D フォールバックなし計算（120胞体 RZ のみ有効）；STL `STLs_v4/` 280枚；vertex-up を 600胞体に適用（打破不可能と結論）；600胞体 w 分布解析 |
| 2026-05-13 | `gen_120cell_faceup_nets.m` v4 fix 適用・全12 C2 対応；`paper_draft.tex` face-up 比較テーブル更新；`unfold4D.m` grayscale・黒エッジ・Neutral lighting；SnubCube mirror 表示追加；Figure 6/7 更新；4D nets figure 追加 |
| 2026-05-14 | 120胞体 RZ nets 幾何学的多様性解析；自己交差構造解析（帯間 71.4% / 帯内 28.6%）；帯境界代替エントリー検証（代替なし）；`paper_draft.tex` 全体レビュー・投稿準備 |
| 2026-05-15 | Discussion 末尾に「Relation to the Companion Implementation」追加；`summary.tex` に対応節追加；コンパニオン論文との実装差異整理 |
| 2026-05-16 | 論文タイトル変更，Introduction 平文化，Section 3.3/3.4/3.5 順序入れ替え，Discussion 副節整理；git 初期化・初期コミット；LICENSE（MIT）・README.md 作成；`STLs_v4/` を git に追加；`paper_draft.tex` は git 追跡から除外（arXiv 公開後に追加予定） |
| 2026-06-09 | Face-rotation BFS ネットの全6正凸4-多胞体検証（773/773 有効）；RegionMember バグ発見・修正；arXiv ドラフト `face_rotation_net_paper.tex` 作成・校正 |
| 2026-06-11 | 3D ランダム多面体実験：Delaunay vs Voronoi 比較（N=30），N スケーリング（N=20–200），Thomson 緩和実験；主結論「成功率は面推移性（対称性）が決定的，幾何的均一性は無関係」 |
| 2026-08-11 | Stella4D 調査（4D ネット表示ソフト，重なり判定は非対応）；`face_rotation_net_paper.tex` に `\bibitem{Webb2026}` と §2.1 の言及を追加；§4 に凸一様64種への拡張を Open Problem として追記（`\bibitem{ConwayGuy1965}`） |
| 2026-08-12 | `uniform4D_wythoff.m` 作成（A₄/B₄/F₄ → H₄）；A₄/B₄/F₄ 30種を全根計算（29種 ALL VALID・x3o4o3x のみ 0/240）；`uniform4D_special.m`（snub 24-cell・grand antiprism）；ridge 検出を O(C²) → 放物型軌道に書き換えて H₄ 到達 |
| 2026-08-13 | `bfs_sat.m` で SAT 移植（19.8→1.5 秒/根）；H₄ 全21,600根計算 → **二分法が反証される**（MIXED 2種）；`uniform4D_prisms.m` で角柱17種；凸一様64種・24,487根が確定；論文を §4 独立節に改稿・タイトル変更・Figure 2 追加；`uniform_nets/` に STL/OFF 67件と `allroots.csv` を公開 |
| 2026-08-15 | `_4DData_uniform/` に64種の組合せデータを公開（5要素形式）；リポジトリのサイズ実測と方針整理 |

### 検討済・対象外

| タスク | 理由 |
|--------|------|
| Face-up 方式の全面切り替え | 600胞体は Impossible のまま，他は改善なし |
| Q1（face-up 成功/失敗 C2 の幾何学的違い） | スコープ外 |
| orbit {4,18,51} の柱状構造の幾何学的理由 | 主結論を変えない |
| 600胞体の構造的ボトルネックへの深い対処 | 多段 BT は非現実的 |
| 2段階選択（RZ→R2）の実装と評価 | 論文の見通しを悪くする |
| 左利きペーリングとの対称性検討 | 結果が鏡像になるだけの可能性が高い |

### 未解決・未実施

| 優先度 | タスク |
|:------:|--------|
| 高 | MIXED 2種（x5o3x3x, x5x3o3x）が**なぜ根に依存するか**の構造分析（2026-08-15 時点で最大の謎） |
| 高 | JCDCG³ 2026 発表内容と改稿後論文の整合（方針：`And more results...` の形で部分的に新情報） |
| 中 | S. Chaidee の ORCID 取得・両論文に追記 |
| 中 | 3つの結果（全根 valid / 根依存 / 全根 invalid）を分ける不変量の探索 |
| 低 | Kaino (2019) との内容面の比較（未実施） |

**2026-08-15 に完了扱いへ移した項目**

| 旧タスク | 現状 |
|---------|------|
| GitHub リポジトリ作成 → push | 完了。<https://github.com/takashi-randomwalker/apple-peel-4d> |
| arXiv に apple-peel ドラフト公開 | 完了。**arXiv:2605.30373** |
| arXiv に face-rotation BFS ドラフト公開 | 投稿済み **submit/7751072**（2026-06-25，cs.CG）。**on hold のまま**。help@ に2回問い合わせ済みだが動かず。JCDCG³ 採択済みのため本筋の損失はないと判断し放置 |

---

## Face-rotation BFS ネット（2026-06-09 セッション）

### 動機

Apple-peel とは独立に，BFS spanning tree + 共有面周り回転という単純な展開法を全6正凸4-多胞体に適用した。Devadoss & Harvey (2022) が 5/8/16-cell で証明した "all spanning trees valid" の open cases（24/120/600-cell）を BFS の範囲で解決することが目的。

### 調査・計算の経緯

1. **120-cell のみ確認**（セッション前半）：`face_rotation_net_120cell.m` で 120/120 有効を確認。BFS layer 0–9 構造も可視化。
2. **600-cell**：`face_rotation_net_600cell.m` を新規作成 → 初回結果 0/600（後でバグ判明）。
3. **5/8/16/24-cell**：`face_rotation_net_small4D.m` を新規作成 → 5/5，8/8 有効，16/16 と 24/24 は 0（後でバグ判明）。
4. **D&H との不整合発見**：D&H が 16-cell "all spanning trees valid" と証明しているのに 0/16 は矛盾 → 詳細調査へ。
5. **バグ発見**（`check_16cell_detail.m`）：
   - RegionMember 判定：62 ペアで "hit"
   - 境界距離検査（Test B）：0 ペア → 全ヒットが境界接触
   - 体積基準（Test C, `RegionMeasure > eps`）：0 ペア
   - **結論**：`RegionMember` は閉領域判定のため境界接触（共有辺・頂点）を誤って重複と扱う
6. **修正版**（`face_rotation_net_all4D_v2.m`）：`trueIntersectQ` を `RegionMeasure[RegionIntersection] > 10^-8` で実装 → 全6胞体 ALL VALID（773/773）。

### バグ影響のまとめ

| 多胞体 | RegionMember（誤） | RegionMeasure（正） |
|--------|:-----------------:|:-------------------:|
| 5-cell   |  5/5  |  5/5 |
| 8-cell   |  8/8  |  8/8 |
| 16-cell  | **0/16** | 16/16 |
| 24-cell  | **0/24** | 24/24 |
| 120-cell | 120/120 | 120/120 |
| 600-cell | **0/600** | 600/600 |

16-cell と 24-cell（正八面体セル，接触多い）および 600-cell（正四面体 600 個，非常に多くの境界接触）が影響を受けた。5-cell・8-cell・120-cell は結果不変。

### 作成スクリプト

| ファイル | 内容 |
|----------|------|
| `face_rotation_net_120cell.m` | 120-cell BFS ネット可視化（TemperatureMap 色付け） |
| `face_rotation_net_600cell.m` | 600-cell 専用 2 フェーズ検証（bbox proxy + 全テスト） |
| `face_rotation_net_small4D.m` | 5/8/16/24-cell 汎用検証 |
| `check_16cell_detail.m` | RegionMember vs 体積基準 比較調査 |
| `face_rotation_net_all4D_v2.m` | 全6胞体 体積基準統合版（**現行推奨**） |
| `face_rotation_net_viewer.nb` | BFS アニメーション付き可視化 Notebook |
| `face_rotation_net_paper.tex` | arXiv ドラフト論文（5 ページ） |

### 論文の内容

タイトル："Face-Rotation BFS Nets of the Six Regular Convex 4-Polytopes"

- Theorem 1（Main）：全6正凸4-多胞体の任意の根から生成した face-rotation BFS ネットは有効
- Algorithm 1（UnfoldTransform）：4D アフィン変換の疑似コード
- 重要な方法論注記：体積基準 vs RegionMember の違い（Definition + Remark）
- Apple-peel 対比表（120-cell: ordering 1440/1440 成功 → net 0/1440，BFS 120/120；600-cell: ordering 0/2400，BFS 600/600）
- Open problems：(1) 理論的証明，(2) 全スパニングツリーへの拡張，(3) 非正則 4-多胞体

### 参考文献（face-rotation 論文）

- Devadoss & Harvey (2022): CGTA, vol 111, article 101977, 2023; arXiv:2111.01359
- Coxeter (1973): Regular Polytopes, Dover
- Shephard (1975): Dürer conjecture
- Turney (1984): Journal of Recreational Mathematics, 17(1):1–16（超立方体 261 ネット）
- Yoshino & Chaidee (2026): arXiv:2605.30373（companion apple-peel paper，2026-06-01 公開済み）

---

## 論文修正履歴

### 2026-05-16 セッション

| 修正 | 詳細 |
|------|------|
| line 671: "four" → "three" | RZ 100% は 3 種（TruncOcta, TruncCubocta, TruncIcosa） |
| line 671-674 文構造整理 | "(all 100%)" 削除，"and is best on" → 完全な文 |
| line 703: "four solids" → "three solids" | Truncated Icosidodecahedron は 66.7% で除外 |
| line 1258-1259 | 16/24-cell の geo-unique=1 確認に基づく根拠追加 |
| line 1386 | "9 symmetric bands" → "9 bands (grouped by symmetry orbit)" |
| line 1394 | `$w \approx 2.118$` 削除（正規化が未定義） |
| Conclusion | SnubCube/TT の partial 結果追加（13 種の内訳 3+1+2+7=13） |
| 4D 条件式 | `\label{eq:4dleft}` 付与（`equation` 環境に変更） |
| line 1094-1095 | `eq:left`（3D）→ `eq:4dleft`（4D），"base point" → "reference plane (c₁,c₂)" |
| Figure 5 キャプション | TT は RZ で表示中だが RS が高い（66.7%）ことを明記 |
| `tab:3dnet` キャプション | 600-cell の em-dash が Impossible を意味することを補足 |
| タイトル | "Apple-Peel Unfolding in Three and Four Dimensions: Spiral and Zonal Selection Rules" |
| Introduction | 箇条書き → 平文散文 |
| Section 2.2 | Remark 2.2 削除；"Algorithm~\ref{alg:peel} gives..." を Remark 2.1 末尾に移動 |
| Section 3.1 Remark 3.1 | 「in this case」「in general cases」修飾語追加；方位角間隔の段落削除 |
| Sections 3.3/3.4/3.5 順序入れ替え | Mirror（新3.3）→ TruncIcosa（新3.4）→ SnubCube（新3.5） |
| Table 3 (tab:mirror) | 非ゼロセルにバリアント数括弧追記（例 `24\,(2)`）；ゼロセル `\multicolumn{1}{c}{-}` |
| Discussion 6.3 削除 | 「120-Cell: Valid Nets and Orientation」を Section 5.4 末に paragraph 移動 |
| Section 4.1 | 参照先を `sec:ex-120cell`（5.4 節）に更新 |
| Section 5.4 末 | `\paragraph{3D-face-centroid-up orientation.}` 追加 |
| Discussion 構成 | 4 節 → 2 節（無題段落 2 つ + 節 2 つ） |
| Conclusion Future work | ランダム凸多面体（凸包の成功率解析）を追加 |

### 2026-05-15 セッション

- Discussion 末尾に「Relation to the Companion Implementation」追加（コンパニオン論文 arXiv:2604.16204 との実装差異 3 点）
- キーワードから "Darboux frame" 削除
- `summary.tex` に「コンパニオン論文との実装差異」節追加（日本語）
- `CLAUDE.md` に同節追加

### 2026-05-14 後半セッション

| 修正 | 詳細 |
|------|------|
| Section 4.2 削除 | "4D Darboux Frame via xy-Projection"（`peeling4Df.m` 旧版の説明） |
| Table 4dcomp 削除 | xy-projection vs global reference 比較表 |
| Observations 内の xy-projection 言及削除 | 24/120/16胞体から旧数値除去 |
| "three-rule" → "two-rule" | 序文 line 139 |
| Dodecahedron "RS" → "any" | tab:summary |
| 5件の未引用文献に引用追加 | Demaine2007（序文），AronovORourke1992（序文），Akitaya2024（序文），Buekenhout1998・Devadoss2022（4D Setup 節末） |
| `\ref{sec:4dresults-global}` 修正 | → `\ref{sec:4dresults}` |
| Dodecahedron 節の誤図参照削除 | fig:peel-example（Truncated Icosahedron 図）を削除 |
| Abstract に ordering vs validity 区別を明記 | 120 胞体が Perfect でも 3D valid net は 0 件 |
| Perfect/Possible/Impossible 定義を前出し | Section 3 冒頭 |
| Discussion 新節「Spiral vs. Zonal Strategy in Four Dimensions」 | 4D で RS が劣化する構造的理由 |
| Kaino (2019) との比較追記 | `\bibitem{Kaino2019}` 追加；各セクション（Dodecahedron, 5-cell, 120-cell, 600-cell）に対比文 |

---

## 3D ランダム多面体実験（2026-06-11 セッション）

### 動機

4D 正多胞体の Delaunay 類（16-cell/600-cell = Possible/Impossible）と Voronoi 類（8-cell/120-cell = Perfect）の差が 3D でも再現するか確認。さらに，反発型確率過程（Thomson 問題）で点の分布を均一化した際の成功率変化を調査。

### 実験 A：Delaunay vs Voronoi 比較（N=30，20 試行）

```mathematica
pts = Table[sphericalRandom[], {N}];
m   = ConvexHullMesh[pts];      (* Delaunay：全三角形，3-正則 *)
p   = DualPolyhedron[m];        (* Voronoi 双対：混合多角形，平均次数≈6-12/N *)
```

| 条件 | RS | RZ |
|------|---:|---:|
| Delaunay（|F|≈2N-4=56） | ~0% | ~1% |
| Voronoi 双対（|F|=N=30） | ~0% | ~30% |

- Voronoi 双対の面次数例（N=30 シード 42）：4〜8，平均 5.6，SD≈1.2

### 実験 B：N スケーリング（Del/Vor × RS/RZ）

`run_random_scaling.m` による結果（nTrials = 30/20/10/5/3）：

| N | Del RS | Del RZ | Vor RS | Vor RZ |
|--:|-------:|-------:|-------:|-------:|
| 20 | ~0% | ~1% | ~0% | ~51% |
| 30 | ~0% | ~1% | ~0% | ~22% |
| 50 | ~0% | ~0% | ~0% | ~9% |
| 100 | ~0% | ~0% | ~0% | ~0.1% |
| 200 | ~0% | ~0% | ~0% | ~0% |

Voronoi RZ のみ N<50 で有意な成功率を示し，N が増えるにつれて単調に 0% へ収束。

### 実験 C：Thomson 緩和（N=30，ステップ数可変）

`run_repulsive_comparison.m`：Coulomb 斥力の勾配降下で球面上の点を均一化。

**実験1（N=30，10 trials/level）：**

| steps | VorRZ | E_C | SD(deg) |
|------:|------:|----:|--------:|
| 0 | 27.9% | 445.6 | 1.21 |
| 10 | 30.6% | 438.5 | 1.22 |
| 50 | 26.2% | 442.4 | 1.16 |
| 200 | 25.5% | 613.1 | 1.28 |
| 500 | 21.0% | 508.0 | 1.27 |
| 2000 | 27.0% | 453.8 | 1.24 |

E_C が非単調なのは実験設計上の問題（各 steps レベルで異なる乱数初期点，10 試行では分散大）。

**実験2（N={20,30,50}，random vs Thomson steps=2000，各 10 試行）：**

| N | random | Thomson |
|--:|-------:|--------:|
| 20 | 44.7% | 43.9% |
| 30 | 32.6% | 22.1% |
| 50 | 5.9% | 6.5% |

差はいずれも統計誤差（σ≈5–10%）の範囲内。

### 主結論

1. **面推移性（対称性）が支配的因子**：Thomson 均一化（点の等間隔化）は成功率を改善しない。正多胞体が Perfect なのは点配置の均一さではなく，面推移群による等変性（0% or 100% 二択）が理由。
2. **スケールが支配的**：N が大きくなるほど成功率は 0% へ収束。面が六角形に近づいても効果なし。
3. **4D との対比**：Delaunay（3-正則）≈ 16-cell，Voronoi（6-正則）≈ 8-cell という対応が 3D でも成立。
4. **論文 Future Work への含意**：「S² ランダム凸包での scaling 解析」の実証的基盤として位置付けられる。

---

## 凸一様4-多胞体への拡張（2026-08-11〜15 セッション）

### 発端

Stella4D（4D 多胞体の 3D ネットを表示できる唯一の主要ソフト）の調査から始まった。マニュアル §15.8 に「in 4D ... intersections between cells are ignored」とあり，**重なり判定は射程外**＝論文の貢献は先取りされていないと確認。引用を追加した際に，§4 Open Problems の「非正則への拡張」を漠然とした一文から**凸一様多胞体64種**という具体的な目標に書き換えた。それが実際にやってみる流れになった。

### 経緯と，覆った2つの見立て

**第1段階（A₄/B₄/F₄ 30種・1,893根）**：29種が全根 valid，runcinated 24-cell (x3o4o3x) のみ全240根 invalid。**MIXED がゼロ**だった。これを「二分法が正則性を超えて現れる」＝最大の発見として記録した。apple-peel 論文の Proposition 3.1（面推移的な立体は 0% か 100%）と同型の現象に見えたため。

**第2段階（H₄ 15種を root 1 のみ）**：15/15 が VALID。全根は 2640セル級で15時間かかると見積もり，「30種すべてで全か無かだったので root 1 は強い指標」と述べた。

**第3段階（H₄ 全21,600根）**：**両方の見立てが誤りだった。**

- **二分法は反証された**。x5o3x3x が 47/2640，x5x3o3x が 2570/2640 で MIXED
- **root 1 は指標にならない**。x5o3x3x は root 1 で VALID だったが，実際には 2640根中47根しか valid でない。root 1 がたまたま当たりだった

x5x3o3x の失敗根は重なりペアが1個ちょうどだったため，SAT の eps（10⁻⁶）が判定を決めている可能性を疑い，`verify_mixed_h4.m` で裏を取った。失敗根5個・成功根3個 × 2多胞体の16根で **SAT と RegionMeasure が判定・重なり数とも完全一致**，eps を 10⁻⁴〜10⁻¹⁰ に振っても重なり数は不変。数値誤差ではなく本物の交差だった。

**教訓**：30種・2,357根という規模でも，そこから外挿した性質（二分法）と近道（root 1 で代表させる）は両方とも外れた。全数計算を実際にやるまでは仮説として扱うべきだった。

### SAT 移植：ボトルネックは予想と違った

H₄ 全根を可能にするため `RegionMeasure[RegionIntersection[...]]` を分離軸判定に置き換えた（`unfold3DExport.m` の実装を移植）。**ところが SAT 化だけでは 19.8 → 10.7 秒/根にしかならず，プロファイルすると SAT 判定自体の所要は 0.04 秒＝実質ゼロだった。**

| 段階 | x5x3x3x（2640セル）1根あたり |
|------|---:|
| RegionMeasure 版 | 19.8 s |
| SAT 化のみ | 10.7 s |
| faceUp / 展開変換の行列積化，BFS visited の Association 化 | 8.5 s |
| 候補ペア抽出を密行列 → 重心の半径検索に置換 | **1.5 s** |

SAT 化で露出した3つのボトルネックを潰して初めて13倍に届いた。とくに候補ペア抽出は C×C の bbox 密行列比較（2640セルで700万要素×18回）が全体の7割を占めていた。「2セルが交わるには重心間距離が R_i + R_j 以下」という条件で `Nearest` の半径検索に置き換えた。

**検証**：`validate_sat.m` で30種1,893根を再計算し，valid 数だけでなく**根ごとの重なりペア数の列**まで一致することを要求。30/30 MATCH。総計 773 → 85 秒（9.1倍）。

### セル次数は答えではなかった

当初，x3o4o3x が唯一の例外である理由を「最小セル次数5，かつ次数5のセルが8割」に求めかけた。600胞体（次数4）が apple-peel で impossible だったことと方向性が一致していたため。しかし

- snub 24-cell・grand antiprism は**最小次数4なのに全根 valid**
- 角柱は蓋セルが**次数92**に達するが全根 valid，一方 MIXED の x5x3o3x も次数32のセルを持つ

で，低次数も高次数も答えではないと判明した。**3つの結果を分ける不変量は未発見。**

### 生成器の設計

| ファイル | 手法 |
|----------|------|
| `uniform4D_wythoff.m` | Coxeter 図 → Gram 行列 → Cholesky で鏡法線；群は反射の閉包；種点は「リング付き鏡から距離1，なしは0」（鏡 i が生む辺長は 2b_i なので等距離＝等辺長）；d-面は階数 d の放物型部分群の軌道 |
| `uniform4D_special.m` | snub 24-cell と grand antiprism を**600胞体からの diminishing** として構成。除いた頂点の隣接頂点のうち生き残った分が新セル（前者は12個＝二十面体，後者は10個＝五角反柱） |
| `uniform4D_prisms.m` | 角柱は組合せ的に直接構成（凸包計算不要）。オイラーは P の標数2から自動的に0 |

**検証は文献値に依存させなかった**。B₄/F₄/H₄ の V/E/F/C を記憶で書き下すのは誤りの元なので，V = \|W\|/\|W_unringed\|，C = Σ \|W\|/\|W_Jᵢ\|，χ = 0，全 facet 超平面が支持面，各 ridge がちょうど2セル，頂点次数の一様性，を群論とオイラーから自前で予測して照合した。48/48 通過。

### 副産物：`_4DData/f*.m` の形式についての発見

- **第3要素は「各頂点の隣接頂点リスト」**であって「セルの頂点リスト」ではない。f5.m では頂点数5・次数4がセル数5・4頂点と偶然一致するため誤読しやすい（実際に一度誤読した）
- **f8.m の第2要素は32本の辺を両方向で64件**列挙している（f5/f16/f24 は1回ずつ）。誰も読まないので実害なし
- **`face_rotation_net_all4D_v2.m` は `raw[[1]], raw[[4]], raw[[6]]` の3つしか読まない**。第3・第5要素（頂点隣接・面隣接）はリポジトリ内のどこからも参照されていない。一方**ネットに必要なセル隣接は形式に含まれておらず**，毎回 `cellsF` から再計算されている

この発見が `_4DData_uniform/` の5要素形式（未使用の2つを削り，セル隣接を足す）につながった。

### 命名の落とし穴

**x5o3x3x と x5x3o3x を一度取り違えた**（論文の草稿で 47 と 2570 を逆に書いた）。セル構成で判別すること。

- x5o3x3x = $t_{0,2,3}\{5,3,3\}$ = runcitruncated **600-cell**：120 truncated icosahedra + 720 pentagonal prisms + 1200 hexagonal prisms + 600 cuboctahedra
- x5x3o3x = $t_{0,1,3}\{5,3,3\}$ = runcitruncated **120-cell**：120 truncated dodecahedra + 720 decagonal prisms + 1200 triangular prisms + 600 cuboctahedra

45種すべてに英語名を付けるのは転記ミスの元なので，公開データの name 列は**リング配置から機械生成した Schläfli t 記法**にした。
