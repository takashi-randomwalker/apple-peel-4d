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
| 高 | GitHub リポジトリ作成 → `git remote add origin` → `git push` |
| 高 | arXiv にドラフト公開（意図的に延期中，2026-05-16） |
| 高 | 査読対応：等変性を Remark → Proposition + Proof（CGTA 投稿前） |
| 高 | Code availability statement を `paper_draft.tex` に追加 |
| 中 | 参考文献の追加：Schlickenrieder, Pak, Bern et al. など |
| 低 | Kaino (2019) との内容面の比較（未実施） |

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
