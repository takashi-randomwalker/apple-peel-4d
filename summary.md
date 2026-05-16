# 4次元正多胞体に対するApple-Peel展開
## アルゴリズム設計と実装のまとめ

**Takashi Yoshino**

---

## 目次

1. [背景と目的](#1-背景と目的)
2. [3次元アルゴリズムの概要](#2-3次元アルゴリズムの概要)
3. [4次元への拡張](#3-4次元への拡張)
4. [左側条件の再検討と4次元への厳密な拡張](#4-左側条件の再検討と4次元への厳密な拡張)
5. [`peeling4Df` の実装](#5-peeling4df-の実装)
6. [ファイル構成と使い方](#6-ファイル構成と使い方)
7. [計算結果](#7-計算結果)
8. [展開図の3次元可視化と検証](#8-展開図の3次元可視化と検証)
9. [1段バックトラック拡張（`peeling1back.m`）](#9-1段バックトラック拡張peeling1backm)
10. [今後の課題](#10-今後の課題)
11. [正多面体における左条件の詳細分析・規則の呼び名](#11-正多面体における左条件の詳細分析2026-05-0607)
- [付録: peeling4Df 完全コード](#付録-peeling4df-完全コード)

---

## 1. 背景と目的

### 1.1 Apple-Peel展開とは

**Apple-Peel展開**（りんごの皮剥き展開）とは，多面体や多胞体の面（セル）を
リンゴの皮を螺旋状に剥くように順番に選択し，展開図を得る手法である．
Yoshino, Chaidee, Sawatdithep によって多面体（3次元）向けに定式化され，
論文 `applePeelingV01.pdf` に記載されている．

本プロジェクトの目的は，この手法を **4次元正多胞体**に拡張することである．

### 1.2 対象とする正多胞体

6種類の正多胞体を対象とする（表1）．

**表1: 正多胞体の基本情報**

| 名称 | セル数 | 面の種類 | 頂点数 | データファイル |
|:---|---:|:---|---:|:---|
| 5胞体 (pentachoron) | 5 | 三角形 | 5 | `f5.m` |
| 8胞体 / 超立方体 (tesseract) | 8 | 正方形 | 16 | `f8.m` |
| 16胞体 (hexadecachoron) | 16 | 三角形 | 8 | `f16.m` |
| 24胞体 (icositetrachoron) | 24 | 正八面体 | 24 | `f24.m` |
| 120胞体 (hecatonicosachoron) | 120 | 正十二面体 | 600 | `f120.m` |
| 600胞体 (hexacosichoron) | 600 | 三角形 | 120 | `f600.m` |

### 1.3 参照文献

- Yoshino, Chaidee, Sawatdithep: *Apple-Peel Unfolding of Polyhedra* (`applePeelingV01.pdf`) — 3D版の定式化
- Kaino (2019): 24胞体，120胞体，600胞体の4D展開図を直交射影を用いて幾何的・手作業で構成した先行研究

---

## 2. 3次元アルゴリズムの概要

### 2.1 Algorithm 1 (applePeelingV01)

**アルゴリズム1: Apple-Peel Unfolding（3次元版）**

```
入力: 多面体の頂点座標 vers，面リスト faces，隣接リスト nbrs，出発面 F1，次面 F2
出力: (F1, F2) から始まる面選択順序 order，成功フラグ success

1: v ← vers − mean(vers)          // 重心を原点に移動
2: ĉ1 ← 面 F1 の重心（単位ベクトル）
3: v ← Rotate(v, ĉ1 → +ẑ)        // F1 を +z 方向へ整列
4: order ← [F1, F2]
5: ĉi ← 全面 i の重心を再計算
6: last ← F2
7: while 未選択の隣接面が存在する do
8:   L ← { j ∈ nbrs[last] | ĉ_last,y · ĉ_j,x − ĉ_last,x · ĉ_j,y ≥ 0 }
                                   // 左半空間候補
9:   if L ≠ ∅ then
10:    next ← argmax_{j∈L} ĉ_j,z  // 左側候補の中で z 最大
11:  else
12:    next ← argmin_{j∈nbrs[last]} ĉ_j,z  // フォールバック：z 最小
13:  end if
14:  order.append(next)
15:  nbrs[i] ← nbrs[i] \ {next}   for all i
16:  last ← next
17: end while
18: return (order, |order| = |faces|)
```

### 2.2 剥きやすさの分類

- **Perfect**: すべての $(F_1, F_2)$ の組み合わせで全面選択に成功
- **Possible**: 一部の $(F_1, F_2)$ で成功，一部で失敗
- **Impossible**: すべての $(F_1, F_2)$ で失敗

Archimedean 立体13種：Perfect 4，Possible 3，Impossible 6 が報告されている．

### 2.3 Mathematica 実装 `peeling3Df`

実装は `221026ArchimedeanSolidsPt02.nb` 中の `peeling3Df` 関数にある．
シグネチャは以下のとおり：

```mathematica
peeling3Df[origVers_, faces_, top_]
```

隣接面の計算には `RotateLeft` による辺ペアのトリックを用いる．
フォールバック（右側候補のみの場合）は，論文の「最小 $z$」ではなく
`neighbors[[lastChosen]][[1]]`（最初の要素）を取る簡略実装になっている．

返り値は各 $F_2$ に対して `{rotatedVers, faceOrder, successBool}` のリスト．

---

## 3. 4次元への拡張

### 3.1 座標系と剥き軸

4次元空間を $O\text{-}xyzw$ と定め，**$w$ 軸**を剥き軸（3Dの $z$ 軸に対応）とする．

### 3.2 4Dデータ構造

データは `_4DData/f[n].m` に格納されており，読み込むと以下の6要素リストが返る：

```
{vers, edgs, neis, faces, fneis, cells}
```

**表2: f[n].m のデータ構造**

| 変数名 | 形状 | 内容 |
|:---|:---|:---|
| `vers` | (頂点数) × 4 | 4次元頂点座標 |
| `edgs` | (辺数) × 2 | 辺を構成する頂点番号のペア |
| `neis` | (頂点数) × * | 各頂点の隣接頂点番号リスト |
| `faces` | (面数) × * | 面を構成する頂点番号リスト |
| `fneis` | (面数) × * | 各面の隣接面番号リスト |
| `cells` | (セル数) × * | **セルを構成する面番号リスト** |

**重要**: `cells[[i]]` は**面番号のリスト**（頂点番号ではない）．
セル $i$ の頂点集合は

$$V_i = \bigcup_{f \in \texttt{cells}[[i]]} \texttt{faces}[[f]]$$

として得られる．

例（5胞体, n=5）：

```mathematica
cells = {{1,2,4,7}, {1,3,5,8}, {2,3,6,9}, {4,5,6,10}, {7,8,9,10}}
(* 各セルは4枚の三角形面からなる四面体 *)
```

### 3.3 3D から 4D へのアナロジー対応

| 3次元 | 4次元 |
|:---|:---|
| 多面体の面 $F$ | 多胞体のセル $C$（3D多面体） |
| 剥き軸 $z$ | 剥き軸 $w$ |
| 面の重心 $\hat{c} \in \mathbb{R}^3$ | セルの重心 $\hat{c} \in \mathbb{R}^4$ |
| $F_1$ を $+z$ 方向へ回転 | $C_1$ を $+w$ 方向へ回転 |
| 左半空間 $(c_k \times \hat{z}) \cdot c' \geq 0$ | 左半空間 $c_y c'_x - c_x c'_y \geq 0$ |
| 隣接 = 辺共有（頂点2点以上共有） | 隣接 = 2D面共有（面番号1つ以上共有） |

### 3.4 左半空間条件の導出

3Dの条件 $(\hat{c}_k \times \hat{z}) \cdot \hat{c}' \geq 0$ を $\hat{z} = (0,0,1)$ で展開すると

$$c_{k,y}\,c'_x - c_{k,x}\,c'_y \geq 0.$$

この式は $c$ と $c'$ の $xy$ 成分のみに依存しており，4次元空間においても**同じ式をそのまま適用できる**：

$$\boxed{c_y\,c'_x - c_x\,c'_y \geq 0}$$

これは，$w$ 軸まわりの回転（$wz$ 平面内での回転）と等価であり，数式的に最も簡潔な拡張形である．

### 3.5 4D回転の実装

$C_1$ の重心 $\hat{c}_1$ を $+w$ 方向に向ける回転には，Mathematicaの `RotationMatrix` を用いる：

```mathematica
angle4D = ArcCos[Clip[cc1 . wHat / Norm[cc1], {-1, 1}]];
RotationMatrix[angle4D, {cc1/Norm[cc1], wHat}] . # & /@ localVers
```

`RotationMatrix[angle, {u, v}]` は任意次元で $u$–$v$ 平面内の回転行列を生成するため，4次元にも直接適用できる．

---

## 4. 左側条件の再検討と4次元への厳密な拡張

本節では，3次元 Apple-Peel における「一番左の胞を選ぶ」という条件を
厳密に数学的に定式化し，4次元への正しい拡張を議論する．

### 4.1 3次元における左側条件の幾何学的意味

3次元多面体の Apple-Peel では，剥き軸（z軸）に平行な方向に胞を整列させた後，
各ステップで「現在の面の重心 $\mathbf{c}_k$ から見て最も左にある隣接面」を選ぶ．
この「左」は $xy$ 平面上の角度として定義される：

$$c_{k,y}\, c'_x - c_{k,x}\, c'_y \;\ge\; 0 \tag{1}$$

ただし $\mathbf{c}' = (c'_x, c'_y)$ は候補面の重心の $xy$ 成分である．

#### Darboux フレームによる解釈

各ステップに局所フレーム（Darboux フレーム）を導入する：

- **法線方向**: $\hat{n}_k = \mathbf{c}_k / |\mathbf{c}_k|$（球面上の位置）
- **前進方向**: $\hat{f}_k = \mathbf{c}_{k-1} \times \mathbf{c}_k \,/\, |\mathbf{c}_{k-1} \times \mathbf{c}_k|$（前のステップから現在への進行方向の「左方向」を与える外積）
- **左方向**: $\hat{L}_k = \hat{n}_k \times \hat{f}_k$

候補面 $\mathbf{c}_j$ が「左半空間」にあるとは $\hat{L}_k \cdot \mathbf{c}_j \ge 0$ のことである．
これは行列式

$$\det\!\bigl(\mathbf{c}_{k-1},\, \mathbf{c}_k,\, \mathbf{c}_j\bigr) \;\ge\; 0 \tag{2}$$

と正確に等価である．実際，

$$\hat{L}_k \cdot \mathbf{c}_j = \frac{\det(\mathbf{c}_{k-1}, \mathbf{c}_k, \mathbf{c}_j)}{|\mathbf{c}_k|\,|\hat{d}_k^\perp|}$$

ただし $\hat{d}_k^\perp$ は $\mathbf{c}_{k-1}$ の $\mathbf{c}_k$ への垂直成分の正規化ベクトルである（スカラー倍は正）．

> **注意**: 3次元では $\mathbf{c}_{k-1} \times \mathbf{c}_k$ が1次元の（一意な）左方向を与える．これが「左が一意に定まる」理由である．

### 4.2 4次元への拡張の困難：直交補空間の次元

$\mathbb{R}^4$ では，2つのベクトル $\mathbf{c}_{k-1}, \mathbf{c}_k$ が張る2次元部分空間の
直交補空間は**2次元**になる（$\mathbb{R}^3$ では1次元）．
したがって外積のような単一の「左ベクトル」は存在せず，「左」の方向は一意に定まらない．

この困難を回避するために，以下の方針をとる：
4次元多胞体の Apple-Peel では，各胞の重心を $xyz$-超平面（$w$ 成分を無視）に射影して
3次元ベクトル $\tilde{\mathbf{c}}_k = (c_{k,x}, c_{k,y}, c_{k,z})$ を得，
この射影空間で行列式条件を適用する．

### 4.3 4次元への正しい拡張：直前ステップを参照した行列式条件

**定義 4.1（4次元 Apple-Peel 左側条件）**

$k \ge 2$ ステップ目において，現在の胞の重心の $xyz$ 射影を $\tilde{\mathbf{c}}_k \in \mathbb{R}^3$，
直前の胞の重心の $xyz$ 射影を $\tilde{\mathbf{c}}_{k-1} \in \mathbb{R}^3$，
候補胞 $j$ の重心の $xyz$ 射影を $\tilde{\mathbf{c}}_j \in \mathbb{R}^3$ とするとき，
胞 $j$ が**左半空間**にあるとは

$$\det\!\bigl(\tilde{\mathbf{c}}_{k-1},\, \tilde{\mathbf{c}}_k,\, \tilde{\mathbf{c}}_j\bigr) \;\ge\; 0 \tag{3}$$

が成立することをいう．

この条件は，$xyz$ 超平面に制限した Darboux フレームの左方向と同値である．

**初期条件の不定性**

$k=1$（最初の選択）では $\tilde{\mathbf{c}}_{k-1}$ がまだ存在しないため，式(3)は適用できない．
現在の実装（`peeling4Df.m`）では，$k=1$ において $\tilde{\mathbf{c}}_{k-1} = \mathbf{0}$ に相当する条件
（2成分行列式 $c_{k,y} c'_x - c_{k,x} c'_y \ge 0$）を使用しており，
これは $xy$ 平面上の $z$ 軸方向ベクトル $(0,0,1)$ に相当する特定の参照ベクトルを仮定している．

### 4.4 三つの選択規則

左半空間条件を満たす複数の候補の中から次の胞を1つ選ぶ規則として，以下の三つが考えられる．
重心の $xyz$ 射影を球面 $S^2$ 上の点と見なすと，各規則の幾何学的意味が明確になる．
また $\hat{n} = \tilde{\mathbf{c}}_k / |\tilde{\mathbf{c}}_k|$ を現在の位置の単位ベクトル，
$\hat{f}$ を Darboux フレームの前進方向，$\varphi_j$ を
$\hat{n}$ から $\tilde{\mathbf{c}}_j$ へ向かう方向と $\hat{f}$ のなす角と定義する．

**定義 4.2（選択規則）**

左半空間条件を満たす候補集合を $\mathcal{L}_k$ とする．

- **規則1（最大角度）**: $j^* = \arg\max_{j \in \mathcal{L}_k} \varphi_j$．
  最も「左回り」の候補を選ぶ．$\varphi_j$ が最大 $= 90°$ に近いほど強い左旋回を意味し，密な螺旋を生む．
- **規則2（最小角度 / loxodrome）**: $j^* = \arg\min_{j \in \mathcal{L}_k} \varphi_j$．
  最も「前方方向に近い左」を選ぶ．これは球面上の **loxodrome（等角螺旋，rhumb line）** の近似に相当する：常に進行方向に対して一定の小さな左傾きを保ちながら前進する．
- **規則3（最大 $w$ / 緯度優先）**: $j^* = \arg\max_{j \in \mathcal{L}_k} \hat{c}_{j,w}$．
  $w$ 座標が最大の候補を選ぶ．これは球面上の「緯度」（極からの距離）を優先的に下降する軌跡に対応し，同じ緯度を一周してから次の緯度へ移行する **緯度別走査（latitude-by-latitude scan）** の連続版と解釈できる．

### 4.5 規則2（loxodrome）と規則3（緯度優先）の幾何学的比較

**3次元球面上での比較**

3次元では，規則2 は球面上の loxodrome（等角螺旋）に対応し，
規則3 は緯度ごとに輪を描いてから次の緯度に降りる軌跡に対応する．

- **loxodrome** は常に子午線に対して一定の角度を保って進む曲線で，球面をどこまでも巻き続けて極点に収束する無限螺旋になる（有限多面体では近似）．
- **緯度優先（規則3）** は各緯度の円を完全に一周するように選択する離散版軌跡で，緯度が固定された区間では東西方向に進み，胞がなくなると次の緯度（小さい $w$）に移る．

**多胞体での差異**

離散的な多胞体では，各ステップで選べる胞の数が限られるため，連続的な意味での loxodrome や緯度周回は近似にすぎない．しかし選択基準の違いにより，成功率が異なる可能性がある．

### 4.6 Kaino (2019) との対応

Kaino (2019) の手法は，球面に外接する正多面体の面を緯度ごとに走査するという構造を持つ．
これは本稿の **規則3（最大 $w$）** に相当する：
最初の胞（top）付近の高い緯度から始め，次第に低い緯度（$w$ が小さい方向）へ進む．

現在の `peeling4Df.m` の実装は，左半空間候補がある場合に
$\hat{c}_{j,w}$ の降順で選択しており，概念的に規則3に対応する．
ただし，厳密な左条件として式(3)ではなく式(1)の2成分版を用いている点が現状の実装との相違である．

### 4.7 2段階選択の提案

**表3: 2段階選択方針（提案）**

| 段階 | 条件 | 内容 |
|:---:|:---|:---|
| 第1段階 | 規則3（最大 $w$）でフィルタ | $\hat{c}_{j,w}$ が最大の候補を試みる（現在の実装に相当） |
| 第2段階 | 規則2（loxodrome）で補完 | 規則3が詰まった場合，loxodrome 条件（最小 $\varphi$）で再選択 |

この方針の根拠：

- 規則3 は定義が簡単で Kaino との整合性がある．
- 規則2 は理論的に興味深い（球面 loxodrome の離散近似）が，多胞体では規則3 の方が直観的に「螺旋状に剥く」動作を再現する．
- 定義はまず簡単な規則を優先し，うまくいかなければより精緻な規則を適用する，という段階的アプローチが実装上も合理的である．

将来的には，規則1・規則2・規則3の成功率を5胞体・8胞体などの小さな多胞体で組織的に比較し，最適な選択規則を同定することが課題である（§10 参照）．

---

## 5. `peeling4Df` の実装

### 5.1 関数シグネチャ

```mathematica
peeling4Df[origVers_, faces_, cells_, top_]
```

引数：

- `origVers`: 頂点の4次元座標リスト（`vers` に対応）
- `faces`: 面の頂点番号リスト
- `cells`: セルの面番号リスト（`f[n].m` の第6要素）
- `top`: 出発セル $C_1$ のインデックス

返り値：各 $C_2$ 候補（$C_1$ の隣接セル）に対して `{rotatedVers, selectionOrder, successBool}` のリスト．

### 5.2 アルゴリズムの流れ

**アルゴリズム2: peeling4Df（4次元版，グリーディ）**

```
入力: 4次元頂点座標 vers，面リスト faces，セルリスト cells，出発セル top
出力: 各 C2 に対して (v', order, success) のリスト

CellVerts(i) := ∪_{f ∈ cells[i]} faces[f]   // セル i の頂点インデックス集合

1:  v ← vers − mean(vers)                    // 重心を原点に移動
2:  ĉ1 ← mean(v[CellVerts(top)])             // 出発セル C1 の重心
3:  v ← Rotate(v, ĉ1 → +ŵ)                  // C1 を +w 方向へ整列
4:  for i = 1 to |cells| do
5:    ĉi ← mean(v[CellVerts(i)])             // 全セルの重心
6:  end for
7:  for i = 1 to |cells| do
8:    adj[i] ← { j ≠ i : cells[i] ∩ cells[j] ≠ ∅ }  // 隣接リスト構築
9:  end for

10: for each C2 ∈ adj[top] do               // 各第2セルについてペーリング
11:   order ← [top, C2]
12:   nbrs ← copy of adj (全エントリから top, C2 を除去)
13:   last ← C2
14:   while nbrs[last] ≠ ∅ do
15:     L ← { j ∈ nbrs[last] | ĉ_last,y · ĉ_j,x − ĉ_last,x · ĉ_j,y ≥ 0 }
16:     if L ≠ ∅ then
17:       next ← argmax_{j∈L} ĉ_j,w         // 左側候補の中で w 座標最大
18:     else
19:       next ← argmin_{j∈nbrs[last]} ĉ_j,w  // フォールバック：右側候補の中で w 最小
20:     end if
21:     order.append(next)
22:     nbrs[i] ← nbrs[i] \ {next}           for all i
23:     last ← next
24:   end while
25:   yield (v, order, |order| = |cells|)
26: end for
```

### 5.3 3D実装との主な相違点

- セル頂点の取得に `Union @@ faces[[ cells[[i]] ]]` が必要（`cells` が面番号リストのため）．
- 隣接判定が頂点共有数 ≥ 3 から面番号共有数 ≥ 1 に変更．
- フォールバックを「右側で $w$ 最小」として**正確に実装**（3D版は `neighbors[[lastChosen]][[1]]` という簡略実装だった）．
- 4D回転に `RotationMatrix[angle, {u, v}]` の4次元版を利用．

---

## 6. ファイル構成と使い方

### 6.1 ファイル一覧

**表4: 主なファイル**

| ファイル | 内容 |
|:---|:---|
| `peeling4Df.m` | アルゴリズム本体（単体で `Get` して使用可） |
| `run4DPeeling.m` | 全6胞体の実行・分類スクリプト |
| `4DPeelingAll.nb` | 上記2ファイルをまとめたNotebook |
| `_4DData/f[n].m` | 各正多胞体のデータ ($n=5,8,16,24,120,600$) |
| `applePeelingV01.pdf` | 3D版アルゴリズムの論文 |

### 6.2 実行方法

**方法 1 — Notebook を使用**

```mathematica
(* 4DPeelingAll.nb を Mathematica で開き，
   Evaluation -> Evaluate Notebook を実行 *)
```

**方法 2 — スクリプトとして実行**

```mathematica
SetDirectory["/path/to/260324Peeling4D"];
Get["peeling4Df.m"];
Get["run4DPeeling.m"];
```

**方法 3 — 個別に呼び出す**

```mathematica
(* 5胞体のデータを読み込む *)
{vers, edgs, neis, faces, fneis, cells} =
    Get["/path/to/_4DData/f5.m"];

(* 出発セル 1 に対してペーリングを実行 *)
result = peeling4Df[N[vers], faces, cells, 1];

(* result[[k]] = {rotatedVers, order, successBool} for k-th F2 *)
result[[1, 3]]   (* True/False *)
result[[1, 2]]   (* cell selection order *)
```

### 6.3 分類の確認方法

```mathematica
(* 全 F1, F2 の結果を取得 *)
ans = Table[peeling4Df[N[vers], faces, cells, top],
            {top, Length[cells]}];

allBools  = Flatten[Map[Last, ans, {2}]];
trueCount = Count[allBools, True];

classification = Which[
  Count[allBools, False] == 0, "Perfect",
  trueCount == 0,              "Impossible",
  True,                        "Possible"
];
```

---

## 7. 計算結果

### 7.1 分類表

`peeling4Df` を全 $(C_1, C_2)$ の組み合わせに対して実行し，
各正多胞体を Perfect / Possible / Impossible に分類した結果を表5に示す．

**表5: 4次元正多胞体の Apple-Peel 展開 分類結果**

| 名称 | セル数 | 隣接数 | $(C_1,C_2)$ 総数 | ユニーク成功数 | 分類 |
|:---|---:|---:|---:|---:|:---|
| 5胞体 (pentachoron) | 5 | 4 | 20 | 20 | **Perfect** |
| 8胞体 (tesseract) | 8 | 6 | 48 | 48 | **Perfect** |
| 16胞体 (hexadecachoron) | 16 | 4 | 64 | 26 | **Possible** |
| 24胞体 (icositetrachoron) | 24 | 8 | 192 | 119 | **Possible** |
| 120胞体 (hecatonicosachoron) | 120 | 12 | 1,440 | 357 | **Possible** |
| 600胞体 (hexacosichoron) | 600 | 4 | 2,400 | — | 計算未完了 |

**「ユニーク成功数」** とは，全 $(C_1, C_2)$ ペアから得られた成功したセル選択順序を重複除去した数である．

### 7.2 3次元結果との比較

3次元版では Archimedean 立体13種について Perfect: 4，Possible: 3，Impossible: 6 と報告されている．
4次元正多胞体では Impossible は現時点で確認されず，6種のうち少なくとも5種が展開可能（Perfect または Possible）である．
これは正多胞体の高い対称性に起因するものと考えられる．

### 7.3 成功率

**表6: $(C_1, C_2)$ ペアに対する成功率（推定下限）**

| 名称 | ユニーク成功率 | 備考 |
|:---|:---|:---|
| 5胞体 | 100% (20/20) | Perfect |
| 8胞体 | 100% (48/48) | Perfect |
| 16胞体 | ≥ 41% (26/64) | Possible |
| 24胞体 | ≥ 62% (119/192) | Possible |
| 120胞体 | ≥ 25% (357/1440) | Possible |
| 600胞体 | — | 未完了 |

Possible の場合，同一の順序が複数の $(C_1, C_2)$ ペアから生成されることもあり得るため，実際の成功ペア数はユニーク成功数以上となる（表中は下限値）．

---

## 8. 展開図の3次元可視化と検証

### 8.1 4D → 3D 展開（`unfold3DExport.m`）

得られたセル選択順序をもとに，各セル（3次元多面体）を3次元空間に連接させた展開図を生成する実装を `unfold3DExport.m` に記述した．手順は以下のとおりである．

1. **SVD 射影**: 各セルの頂点群（4次元）を特異値分解で最良3次元超平面に等長射影し，局所3次元座標を得る．
2. **Procrustes 整列**: 隣接するセル間の共有面の頂点座標が一致するよう，回転行列 $R$ と平行移動 $\boldsymbol{t}$ を最適化する（スケールなし）．
3. **鏡映補正**: Procrustes 整列後，前セルと新セルが共有面に対して同じ側にある場合，共有面を通した鏡映によって「外側」へ展開する．
4. **STL 出力**: 展開されたセル群を三角形分割し，`Graphics3D` + `Export` で STL ファイルとして保存する．

### 8.2 重なり判定（Separating Axis Theorem）

3次元に展開した際に異なるセル同士が重なるかどうかを，分離軸定理（SAT）を用いて全ペアに対して検査する．

- 候補軸：各セルの面法線ベクトルおよび辺ベクトルの外積
- AABB（軸並行バウンディングボックス）による事前フィルタで高速化
- `eps` $= 10^{-6}$ により共有面での接触は重なりと見なさない

**表7: STL 出力数（`STLs/` フォルダ）**

| 多胞体 | 出力 STL 数 |
|:---|---:|
| 5胞体 | 20 |
| 8胞体 | 48 |
| 16胞体 | 26 |
| 24胞体 | 119 |
| 120胞体 | 357 |
| 600胞体 | — (未完了) |
| **合計** | **570** |

STL ファイルには重なりのある展開図も含まれる．重なりなし（valid net）の数については今後の確認事項である（§10 参照）．

---

## 9. 1段バックトラック拡張（`peeling1back.m`）

グリーディ・アルゴリズムでは，左側候補もフォールバック候補もない状態（行き詰まり）が生じた場合に展開が失敗する．
この問題に対処するため，**1段バックトラック**拡張を実装した（`peeling1back.m`）．

### 9.1 アルゴリズム

グリーディ版と同じ前処理（重心移動・回転・重心計算・隣接リスト構築）を行った後，以下の選択ループを実行する．

**アルゴリズム3: peeling1back（1段バックトラック版，選択ループ部分）**

```
入力: セル重心配列 ĉ，隣接リスト adj，出発セル top，第2セル C2
出力: (order, success)

1:  order ← [top, C2]
2:  nbrs ← copy of adj（top, C2 を全エントリから除去）
3:  last ← C2
4:  justBT ← False
5:  tried ← {}               // 位置 → 試行済みセルの集合
6:  while |order| < |cells| do
7:    pos ← |order|
8:    cands ← nbrs[last] \ tried[pos]   // 未試行の候補
9:    if cands = ∅ then                 // 行き詰まり
10:     if justBT or |order| ≤ 2 then
11:       break                         // 2段以上のバックトラックは行わない
12:     end if
13:     bad ← order.last()
14:     order.removeLast()
15:     nbrs ← 直前の前進前の状態を復元
16:     tried[|order|] ← tried[|order|] ∪ {bad}
17:     justBT ← True
18:   else                              // 前進
19:     L ← { j ∈ cands | ĉ_last,y · ĉ_j,x − ĉ_last,x · ĉ_j,y ≥ 0 }
20:     if L ≠ ∅ then
21:       next ← argmax_{j∈L} ĉ_j,w
22:     else
23:       next ← argmin_{j∈cands} ĉ_j,w
24:     end if
25:     nbrs の現在状態を保存
26:     order.append(next)
27:     nbrs[i] ← nbrs[i] \ {next}      for all i
28:     last ← next
29:     justBT ← False
30:   end if
31: end while
32: return (order, |order| = |cells|)
```

この拡張により，グリーディ版で失敗していた一部の $(C_1, C_2)$ ペアが成功する可能性がある（詳細な比較は今後の課題）．

---

## 10. 今後の課題

1. **600胞体の計算完了**
   600胞体は計算量が大きく（600セル，$(C_1,C_2)$ 総数 2,400），現時点で完全な分類が得られていない．分割実行あるいは計算資源の増強が必要である．

2. **重なり判定の統計整理**
   570枚の STL について，重なりなし（valid net）の比率を多胞体ごとに集計し，表として整理する．

3. **1段バックトラック拡張との比較**
   `peeling1back.m` による成功率とグリーディ版 `peeling4Df` との比較を定量化する．

4. **Kaino (2019) との比較**
   本アルゴリズムが生成した展開図が，Kaino の先行研究（直交射影による手作業構成）と整合するかを確認する．

5. **左利きペーリングの対称性検討**
   $c_y c'_x - c_x c'_y \leq 0$ による左利きペーリングとの結果の対称性・双対性を確認する．

6. **選択規則の体系的比較**
   §4 で提案した規則1（最大角度），規則2（loxodrome），規則3（最大 $w$）の三つの選択規則を，5胞体・8胞体など小規模な正多胞体に適用して成功率を比較する．3次元多面体（Archimedean 立体など）での検証も有効である．

7. **厳密な4次元左条件の実装**
   現在の `peeling4Df.m` が用いる2成分左条件 $c_{k,y}c'_x - c_{k,x}c'_y \ge 0$ を，定義4.1の行列式条件 $\det(\tilde{\mathbf{c}}_{k-1}, \tilde{\mathbf{c}}_k, \tilde{\mathbf{c}}_j) \ge 0$ に置き換えて成功率の変化を確認する（初期条件 $k=1$ の扱いに注意）．

8. **2段階選択の実装と評価**
   表3の2段階選択方針（規則3 優先 → 規則2 補完）を実装し，グリーディ版および1段バックトラック版と成功率を比較する．

---

## 付録: peeling4Df 完全コード

```mathematica
peeling4Df[origVers_, faces_, cells_, top_] :=
  Module[
    {wHat, cc1, localVers, cCenters,
     neighbors, neighborsInit,
     order, cansInit, chosenFrom1,
     lastChosen, newChoice, cans,
     angle4D, rCans, cellVerts},

    wHat = {0, 0, 0, 1};

    (* セル i の頂点番号を返す補助関数 *)
    cellVerts[i_] := Union @@ faces[[ cells[[i]] ]];

    (* Step 1: 重心を原点に移動 *)
    localVers = (# - Mean[origVers]) & /@ origVers;

    (* Step 2: 出発セル C1 の重心 *)
    cc1 = Mean[ localVers[[ cellVerts[top] ]] ];

    (* Step 3: +w 方向へ回転 *)
    localVers = Which[
      Abs[(cc1/Norm[cc1]).wHat - 1] < 0.0001, localVers,
      Abs[(cc1/Norm[cc1]).wHat + 1] < 0.0001, -localVers,
      True,
        angle4D = ArcCos[Clip[cc1.wHat/Norm[cc1], {-1,1}]];
        RotationMatrix[angle4D, {cc1/Norm[cc1], wHat}] . # & /@ localVers
    ];

    (* Step 4: 全セル重心 *)
    cCenters = Mean[localVers[[ cellVerts[#] ]]] & /@ Range[Length[cells]];

    (* Step 5: 隣接判定（面番号共有 >= 1） *)
    neighbors = Table[
      Select[Range[Length[cells]],
        Function[j, j!=i &&
          Length[Intersection[cells[[i]], cells[[j]]]] >= 1]],
      {i, Length[cells]}];

    (* Step 6: C1 を隣接リストから除外 *)
    neighbors = neighborsInit = Complement[#, {top}] & /@ neighbors;

    (* Step 7: C2 の候補 *)
    cansInit = neighbors[[top]];

    (* Step 8: 各 C2 に対してペーリング *)
    Table[
      order = {top}; neighbors = neighborsInit;
      chosenFrom1 = cansInit[[i]];
      AppendTo[order, chosenFrom1];
      neighbors = Complement[#, {chosenFrom1}] & /@ neighbors;
      lastChosen = chosenFrom1; newChoice = lastChosen;

      While[
        Length[Flatten[neighbors]] != 0 &&
        Length[neighbors[[newChoice]]] != 0,

        (* 左半空間: c_y * c'_x - c_x * c'_y >= 0 *)
        cans = Select[neighbors[[lastChosen]],
          Function[j, With[{c = cCenters[[lastChosen]]},
            c[[2]]*cCenters[[j]][[1]] -
            c[[1]]*cCenters[[j]][[2]] >= 0]]];

        newChoice = Which[
          Length[cans] > 1,
            cans[[ First[First[
              Position[cCenters[[cans,4]], Max[cCenters[[cans,4]]]]
            ]] ]],
          Length[cans] == 1, cans[[1]],
          True,  (* フォールバック: 右側で w 最小 *)
            With[{rc = neighbors[[lastChosen]]},
              rc[[ First[First[
                Position[cCenters[[rc,4]], Min[cCenters[[rc,4]]]]
              ]] ]]]
        ];

        lastChosen = newChoice;
        AppendTo[order, newChoice];
        neighbors = Complement[#, {newChoice}] & /@ neighbors;
      ];

      {localVers, order, Length[order] == Length[cells]}
      , {i, Length[cansInit]}]
  ]
```

---

## 11. 正多面体における左条件の詳細分析（2026-05-06〜07）

### 11.1 左条件の参照点：局所 vs 大域

#### 旧実装（局所参照）

旧実装では左半空間条件の参照点として直前の面重心 $\mathbf{c}_{k-1}$ を用いていた：

$$\det(\mathbf{c}_{k-1},\, \mathbf{c}_k,\, \mathbf{c}_j) \;\ge\; 0$$

この条件は「局所的には左」を意味するが，参照点がステップごとに変わるため，
螺旋の方向が各ステップで変動し，大域的な一貫性が保証されない．

#### 新実装（大域参照）

左半空間条件の参照点を開始面 $\mathbf{c}_1$（北極，Face-up 後の固定点）に変更した：

$$\det(\mathbf{c}_1,\, \mathbf{c}_k,\, \mathbf{c}_j) \;\ge\; -\varepsilon, \qquad \varepsilon = 10^{-10}$$

$\mathbf{c}_1$ は全ステップで固定されるため，$z$ 軸周りの螺旋方向が大域的に保証される．

#### 単一候補優先の廃止

旧実装には「候補が1つのときは左条件を適用せずそのまま選択する」単一候補優先があったが，
大域参照と矛盾するため廃止した．全候補数によらず左条件を適用する．

### 11.2 等変性（equivariance）の理論と検証

正多面体は**面推移的**（face-transitive）であり，回転対称群 $G$ が任意の $(F_1, F_2)$ ペアを別のペアに写す．
アルゴリズムが $G$ の下で等変であれば，全ペアの成否が一致し，成功率は 0% か 100% のいずれかになる．

等変性の証明：対称 $\sigma \in G$ に対して合成回転 $A = \mathrm{FaceUp}(\sigma F_1) \circ \sigma \circ \mathrm{FaceUp}(F_1)^{-1}$ を定義すると，変換後の重心は $\mathbf{c}'_j = A\mathbf{c}_j$ となる．$A$ は固有回転（$\det A = 1$）なので

$$\det(\mathbf{c}'_1, \mathbf{c}'_k, \mathbf{c}'_j) = \det(A) \cdot \det(\mathbf{c}_1, \mathbf{c}_k, \mathbf{c}_j) = \det(\mathbf{c}_1, \mathbf{c}_k, \mathbf{c}_j)$$

左候補集合は $\sigma$ によるラベル置換の下で不変であり，Darboux 角 $\varphi$ や $z$ 座標も内積量として $A$ で保存される．したがって R1・R2・R3 の全規則が等変である．

**数値的注意事項**：等変性を数値的に保証するには以下の2点が必要：

1. **eps 閾値 $\varepsilon = 10^{-10}$**：正十二面体の黄金比座標では `N[...]` による有限精度化で $\det \approx 0$ が $\sim 10^{-16}$ の誤差を持ち，符号が対称等価なペア間で逆転する．$\varepsilon$ でこれを吸収する．
2. **R3 の二次タイブレーク**：$z$ 値が $10^{-10}$ 以内で同値な候補は $\arg\max\varphi$（R1）で二次選択する．面インデックス順での単純タイブレークは面ラベルの恣意性に依存し，等変性を破る．

**検証結果**（2026-05-06 実行）：

| 多面体 | Pairs | R1(w) | R1(n) | R2(w) | R2(n) | R3(w) | R3(n) |
|--------|------:|------:|------:|------:|------:|------:|------:|
| Tetrahedron | 12 | 100% | 100% | 100% | 100% | 100% | 100% |
| Cube | 24 | 100% | 100% | 100% | 0% | 100% | 100% |
| Octahedron | 24 | 100% | 100% | 100% | 0% | 100% | 100% |
| Dodecahedron | 60 | 100% | 100% | 100% | 0% | 100% | 100% |
| Icosahedron | 60 | 100% | 100% | 0% | 0% | 100% | 100% |

w = with fallback，n = no fallback．全エントリが 0% または 100% — 等変性理論と一致．

### 11.3 eps の符号による成功率の変化

左条件 $\det(\mathbf{c}_1, \mathbf{c}_k, \mathbf{c}_j) \ge \text{thresh}$ において閾値の符号を変えて比較した（$\varepsilon = 10^{-10}$）：

- **eps > 0（thresh = $-\varepsilon$，現在の設定）**：$\det \approx 0$ の面を左候補に含める（許容あり）
- **eps = 0（thresh = $0$）**：厳密ゼロ境界
- **eps < 0（thresh = $+\varepsilon$）**：$\det \approx 0$ の面を左候補から除外（より厳しい条件）

**成功ペア数の比較（スクリプト `run_platonic_eps_comparison.m` による）**：

**R1 (max φ)**

| 多面体 | Pairs | w: eps>0 | w: eps=0 | w: eps<0 | n: eps>0 | n: eps=0 | n: eps<0 |
|--------|------:|---------:|---------:|---------:|---------:|---------:|---------:|
| Tetrahedron | 12 | 12 | 12 | 12 | 12 | 12 | 12 |
| Cube | 24 | 24 | 24 | 24 | 24 | 22 | 0 |
| Octahedron | 24 | 24 | 24 | 24 | 24 | 18 | 0 |
| Dodecahedron | 60 | 60 | 60 | 60 | 60 | 34 | 0 |
| Icosahedron | 60 | 60 | 60 | 60 | 60 | 35 | 0 |

**R2 (loxodrome)**

| 多面体 | Pairs | w: eps>0 | w: eps=0 | w: eps<0 | n: eps>0 | n: eps=0 | n: eps<0 |
|--------|------:|---------:|---------:|---------:|---------:|---------:|---------:|
| Tetrahedron | 12 | 12 | 12 | 12 | 12 | 12 | 12 |
| Cube | 24 | 24 | 22 | 24 | 0 | 0 | 0 |
| Octahedron | 24 | 24 | 22 | 24 | 0 | 2 | 0 |
| Dodecahedron | 60 | 60 | 48 | 0 | 0 | 0 | 0 |
| Icosahedron | 60 | 0 | 19 | 0 | 0 | 0 | 0 |

**R3 (max z)**

| 多面体 | Pairs | w: eps>0 | w: eps=0 | w: eps<0 | n: eps>0 | n: eps=0 | n: eps<0 |
|--------|------:|---------:|---------:|---------:|---------:|---------:|---------:|
| Tetrahedron | 12 | 12 | 12 | 12 | 12 | 12 | 12 |
| Cube | 24 | 24 | 24 | 24 | 24 | 22 | 0 |
| Octahedron | 24 | 24 | 24 | 24 | 24 | 18 | 0 |
| Dodecahedron | 60 | 60 | 60 | 60 | 60 | 34 | 0 |
| Icosahedron | 60 | 60 | 60 | 60 | 60 | 35 | 0 |

#### 考察

**eps > 0（現在の設定）が唯一の正解**：

- eps > 0：全エントリが 0% または 100% → 等変性成立 ✓，成功率も最大
- eps = 0：中途半端な値（22/24，48/60 など）→ 等変性が破れる ✗．$\det \approx 0$ の面の数値符号が対称等価なペア間でランダムに逆転するため
- eps < 0：全エントリが 0% または 100% → 等変性は成立 ✓ だが，with fallback でも R2 の一部が 0% に落ち，no-fallback は全て 0% になる

#### eps < 0 の幾何学的解釈

thresh = $+\varepsilon$ のとき，$\det(\mathbf{c}_1, \mathbf{c}_k, \mathbf{c}_j) \approx 0$ の面は左候補から除外される．

$\det(\mathbf{c}_1, \mathbf{c}_k, \mathbf{c}_j) = 0$ は $\mathbf{c}_1$（北極），$\mathbf{c}_k$（現在の面），$\mathbf{c}_j$（候補面）が同一経線面上に並ぶことを意味する．すなわち「螺旋の進行方向そのまま，南極方向へ真直ぐ降りる面」が除外される．

- これらの面は fallback に回る
- **R3 の fallback**（$\arg\min z$）はこの「真直ぐ南に降りる面」を最低 $z$ として選ぶ → 同じ面が別ルートで選ばれ，with fallback では eps > 0 と同じ結果になる
- **R1 の fallback**（$\arg\max |\varphi|$）はこの面を選ばない（$\varphi \approx 0$ なので優先度最低）→ 別の面を選ぶが，R1 は元々 $\varphi \approx 0$ の面を好まないため with fallback では結果に影響しない
- **no-fallback** では除外された面が単一候補のとき即終了するため，全て 0% になる

### 11.4 R1 と R3 の等価性（正多面体）

#### 観察

上記の表から，R1 と R3 は：
- with fallback の全ケースで**全く同じ成功率**
- no-fallback の全ケースでも**全く同じ成功率**
- eps の符号によらず常に一致

#### 幾何学的理由：緯度帯構造

Face-up 回転後，正多面体の面重心は「緯度帯（latitude band）」に分かれて等間隔に配置される：

```
Dodecahedron の例:
  北極  {1面}   ← 緯度帯1
  上帯  {5面}   ← 緯度帯2（5面が等間隔）
  下帯  {5面}   ← 緯度帯3
  南極  {1面}   ← 緯度帯4
```

同一緯度帯内の面は $z$ 座標が等しく，隣接面との位置関係が周期的に対称である．

- **R1（max φ）**：最も強く左に曲がる面 = 同一緯度帯で次の面（周回方向）
- **R3（max z）**：最も高い $z$ の面 = 同一緯度帯の面（同じ高さで次の面）

したがって同一緯度帯を周回する間は R1 と R3 が**同じ面を選ぶ**．帯の切り替えも同様のロジックで一致する．

#### Fallback の違い

Fallback は選択規則が唯一異なる場面：

| 規則 | Fallback | 対象 |
|------|----------|------|
| R1 | $\arg\max \|\varphi\|$ | 最も大きく曲がる面（左右問わず） |
| R3 | $\arg\min z$ | 最も低緯度の面 |
| R2 | $\arg\min \|\varphi\|$ | 最も直進に近い面 |

R3 の fallback（min z）は通常選択（max z）と**逆方向**：「現帯で詰まったら最南端へジャンプ」という補完的構造．R1 の fallback は通常選択と一貫（角度最大）．

#### 等価性が崩れるケース

正多面体では等価だが，複雑な多面体では割れる：

| 対象 | R1 | R3 | 差の原因 |
|------|----|----|---------|
| 切頂二十面体（soccer ball） | 80.6% | 35.6% | 五角形・六角形の不規則配置で $\varphi$ と $z$ が一致しない |
| 16胞体（4D） | 18.8% | 40.6% | 4次元の緯度帯が不等間隔で $w$ と $\varphi^{xy}$ が乖離 |
| 120胞体（4D） | ≈0.1% | 24.8% | 大規模・高対称での差が顕著 |

正多面体の高い対称性（等間隔な緯度帯構造）が R1 と R3 を偶然一致させており，より複雑な図形では R3 が優位になる．

### 11.6 eps < 0（厳密条件）の正式計算結果

前節 11.3 では `run_platonic_eps_comparison.m` によって eps < 0 の成功数を確認した．
本節ではこれを RS・RZ のみに絞った正式な計算スクリプト `run_platonic_strict.m` を実行し，
結果を `dataPlatonic_strict.mx`（変数 `platonicDataStrict`）に保存した（2026-05-07）．

#### 設定

- thresh = $+\varepsilon = +10^{-10}$（厳密左条件: $\det(\mathbf{c}_1, \mathbf{c}_k, \mathbf{c}_j) \ge +10^{-10}$）
- R2 を除く RS（max φ）と RZ（max z）の2規則のみ
- with fallback / no fallback の両版を計算

#### 計算結果

| 多面体 | Pairs | RS(w) | RS(n) | RZ(w) | RZ(n) |
|--------|------:|------:|------:|------:|------:|
| Tetrahedron | 12 | 12 | 12 | 12 | 12 |
| Cube | 24 | 24 | 0 | 24 | 0 |
| Octahedron | 24 | 24 | 0 | 24 | 0 |
| Dodecahedron | 60 | 60 | 0 | 60 | 0 |
| Icosahedron | 60 | 60 | 0 | 60 | 0 |

w = with fallback，n = no fallback

#### 考察

1. **等変性は維持**：全エントリが 0% または 100% → 厳密条件でも等変性成立 ✓
2. **with fallback = 100%（Tetra 以外も）**：
   - $\det \approx 0$ の面（経線面上に並ぶ「真南面」）が左候補から外れ fallback に回る
   - RS fallback（$\arg\max |\varphi|$）と RZ fallback（$\arg\min z$）いずれもこの面を別ルートで選ぶため，最終的に eps > 0 と同じペーリング順序が得られる
3. **no fallback = 0%（Tetra 以外）**：
   - 正四面体は「真南面」が存在しないため fallback が不要 → no fallback でも 100%
   - 他の正多面体では「真南面」がある段階で左候補が空になり即終了
4. **RS と RZ が完全一致**：正多面体の緯度帯構造（Section 11.4）が厳密条件でも維持される

#### 参照先（データ構造）

```mathematica
Get["dataPlatonic_strict.mx"]   (* platonicDataStrict を読み込む *)
platonicDataStrict["Dodecahedron"]["RS"]["withFallback"]
(* -> {{top, f2, order, success}, ...} *)
```

トップレベル構造：
```
platonicDataStrict[name] = <|
  "nFaces"      -> Integer,
  "totalPairs"  -> Integer,
  "RS" -> <| "withFallback" -> {...}, "noFallback" -> {...} |>,
  "RZ" -> <| "withFallback" -> {...}, "noFallback" -> {...} |>
|>
```

### 11.7 RS・RZ における Fallback の役割：南極面問題

#### 南極面の幾何学的性質

Face-up 回転後，開始面 F1 の重心 $\mathbf{c}_1$ は $+z$ 軸上にある（$\mathbf{c}_1 = (0,0,z_1)$）．
F1 に対して「反対側」の面（**南極面**）の重心 $\mathbf{c}_\text{south}$ も $-z$ 軸上にある（$\mathbf{c}_\text{south} = (0,0,-z_2)$）．

$\mathbf{c}_1$ と $\mathbf{c}_\text{south}$ はともに $z$ 軸上の点であるため，2ベクトルとして線形従属（反平行）である．
スカラー三重積（行列式）は，線形従属なベクトルを含むと恒等的にゼロになる：

$$\det(\mathbf{c}_1,\, \mathbf{c}_k,\, \mathbf{c}_\text{south})
= \mathbf{c}_1 \cdot (\mathbf{c}_k \times \mathbf{c}_\text{south}) = 0
\quad \text{（$\mathbf{c}_k$ の位置によらず）}$$

これは $\mathbf{c}_k \times \mathbf{c}_\text{south}$ の方向が常に $xy$ 平面内に収まり，$z$ 軸上の $\mathbf{c}_1$ と直交するためである（具体的には $\mathbf{c}_k = (x,y,z)$ として $\mathbf{c}_k \times (0,0,-z_2) = (-yz_2,\, xz_2,\, 0)$ となり，$(0,0,z_1)$ との内積は 0）．

**結論**：南極面は任意のステップ $k$，任意の $\mathbf{c}_k$ に対して $\det = 0$ が成立し，
厳密条件（eps < 0，thresh = $+\varepsilon$）のもとでは**全ステップで左候補から除外される**．

#### no-fallback の失敗は常に「最後の1面」

`dataPlatonic_strict.mx`（RS/RZ no-fallback）の失敗ペアで実際に詰まったステップを確認した結果：

| 多面体 | 面数 $n_F$ | 詰まった order 長さ | 残り面数 |
|--------|:----------:|:------------------:|:-------:|
| Cube | 6 | 5 | **1** |
| Octahedron | 8 | 7 | **1** |
| Dodecahedron | 12 | 11 | **1** |
| Icosahedron | 20 | 19 | **1** |

すべて残り1面（= 南極面のみ）の段階で詰まっている．
南極面は全ステップで左候補から除外されるが，**他の面が存在する間は代替候補が必ずあるため**中間では詰まらない．
南極面だけが残った最後のステップで初めて左候補が空になり，no-fallback は即終了する．

#### Tetrahedron のみ no-fallback でも 100% の理由

正四面体は Face-up 後，残り3面の重心が**すべて同一緯度帯**に等間隔（120° おき）に並ぶ．
いずれの面も $z$ 軸上にはなく，かつ隣接する2面は互いに 120° 離れているため，任意のステップ $k$ で

$$\det(\mathbf{c}_1,\, \mathbf{c}_k,\, \mathbf{c}_j) \neq 0 \quad (\text{全候補 } j)$$

が成立する（南極面に相当する面が存在しない）．
したがって厳密条件でも左候補が空になることがなく，no-fallback でも 100% を達成する．

#### Fallback が担う役割

with-fallback が 100% を保てる理由は，fallback が南極面を「別ルートで」選ぶためである：

| 規則 | Fallback の選択基準 | 南極面を選ぶか |
|------|---------------------|:-------------:|
| RS | $\arg\max \|\varphi\|$（最大絶対角） | **選ばない**（$\varphi \approx 0$ のため最低優先） |
| RZ | $\arg\min z$（最低緯度） | **選ぶ**（南極面は最低 $z$ ） |

実は RS fallback は南極面を**選ばない**．南極面は $\varphi \approx 0$（真正面，左にも右にも曲がらない）なので $|\varphi|$ は最小となり，優先度が最低である．
一方 RZ fallback（$\arg\min z$）は南極面を最低 $z$ として選ぶ．

では RS の with-fallback はなぜ 100% か？
南極面が fallback によって「選ばれないまま」でも，RS は南極面以外の候補から別の面を選び続け，最終的に南極面だけが残ったとき初めて fallback が南極面を pending の唯一の候補として選ぶ（このとき pending = {南極面} のみで fallback に入るため，選択の余地がなく南極面が選ばれる）．

**まとめ：Fallback の本質的役割は「南極面を救済すること」ではなく，左候補が枯渇したときに停止を回避することであり，結果として南極面が最終的に拾われる．**

#### R2 の no-fallback 失敗との対比

RS/RZ と異なり，R2（loxodrome，標準条件）の no-fallback 失敗は中盤で起きる：

| 多面体 | 面数 $n_F$ | 詰まった order 長さ | 残り面数 |
|--------|:----------:|:------------------:|:-------:|
| Cube | 6 | 4 | **2** |
| Octahedron | 8 | 5 | **3** |
| Dodecahedron | 12 | 5 | **7** |
| Icosahedron | 20 | 9 | **11** |

R2 は「最小角（直進方向）」を選び続けるため，pending の全候補が右半空間（$\det < -\varepsilon$）に入る状況を作り出す．
南極面の det = 0 とは無関係に，中盤で右方向にしか候補がなくなって詰まる．
これは RS・RZ の「南極面問題」とは本質的に異なるメカニズムである．

### 11.5 規則の呼び名について

今後 R2（loxodrome）は検討対象から外し，R1 と R3 の2規則を主に扱う．
呼び名の使い分けは以下の通り：

**R2 を含む3規則の比較文脈**では従来通り R1 / R2 / R3 を使う．

**R2 を除いた2規則のみを扱う文脈**では以下の呼び名を使う：

| 記号 | 正式名 | 略称 | 選択基準 |
|------|--------|------|---------|
| R1 | **Spiral rule** | **RS** | max φ（最大左旋回）→ 密な螺旋軌跡 |
| R3 | **Zonal rule** | **RZ** | max z（最高緯度優先）→ 緯度帯を順に走査 |

命名の根拠：
- **Spiral**：最大角度選択が z 軸周りの密な螺旋（coil）を生む
- **Zonal**：同一緯度帯（zone）を一周してから次帯へ降りる帯状走査．気象学・地理学で「緯度帯」を意味する標準用語

---

## 12. アルキメデス多面体の再計算結果（2026-05-07）

### 12.1 背景と設定

旧版（`archimedean_faceup_results.mx`）はアルゴリズム改訂前（局所参照・単一候補優先あり）の結果であり，正多面体の等変性検証と矛盾する実装で計算されていた．
改訂版 `peeling3DLoxo.m`（大域 c_1 参照・eps = $10^{-10}$・単一候補優先なし）で `run_archimedean_faceup.m` を再実行し，結果を `archimedean_faceup_results.mx` に上書き保存した（2026-05-07）．

### 12.2 計算結果

| 多面体 | F | Pairs | RS | RS% | R2 | R2% | RZ | RZ% |
|--------|:-:|------:|---:|----:|---:|----:|---:|----:|
| TruncatedTetrahedron | 8 | 36 | 12 | 33.3 | 12 | 33.3 | 12 | 33.3 |
| Cuboctahedron | 14 | 48 | 0 | 0 | 0 | 0 | 0 | 0 |
| TruncatedCube | 14 | 72 | 0 | 0 | 0 | 0 | 0 | 0 |
| TruncatedOctahedron | 14 | 72 | 72 | **100** | 48 | 66.7 | 72 | **100** |
| Rhombicuboctahedron | 26 | 96 | 0 | 0 | 0 | 0 | 0 | 0 |
| TruncatedCuboctahedron | 26 | 144 | 144 | **100** | 72 | 50.0 | 144 | **100** |
| SnubCube | 38 | 120 | 44 | 36.7 | 0 | 0 | 44 | 36.7 |
| Icosidodecahedron | 32 | 120 | 0 | 0 | 0 | 0 | 0 | 0 |
| TruncatedDodecahedron | 32 | 180 | 0 | 0 | 0 | 0 | 0 | 0 |
| TruncatedIcosahedron | 32 | 180 | 180 | **100** | 60 | 33.3 | 180 | **100** |
| Rhombicosidodecahedron | 62 | 240 | 0 | 0 | 0 | 0 | 0 | 0 |
| TruncatedIcosidodecahedron | 62 | 360 | 360 | **100** | 0 | 0 | 240 | 66.7 |
| SnubDodecahedron | 92 | 300 | 0 | 0 | 0 | 0 | 0 | 0 |

### 12.3 旧結果との比較（主要な変化）

| 多面体 | RS(旧) | RS(新) | R2(旧) | R2(新) | RZ(旧) | RZ(新) |
|--------|-------:|-------:|-------:|-------:|-------:|-------:|
| TruncatedOctahedron | 75.0% | **100%** | 6.9% | 66.7% | 29.2% | **100%** |
| TruncatedCuboctahedron | 0% | **100%** | 1.4% | 50.0% | 0% | **100%** |
| SnubCube | 16.7% | 36.7% | **35.0%** | **0%** | 2.5% | 36.7% |
| TruncatedIcosahedron | 80.6% | **100%** | 0% | 33.3% | 35.6% | **100%** |
| TruncatedIcosidodecahedron | 0% | **100%** | 0% | 0% | 0% | 66.7% |

### 12.4 主な知見

#### RS と RZ の完全一致（12/13種）

正多面体と同様，アルキメデス多面体でも RS と RZ は 12/13 種で同一の成功率を示した．
唯一の例外は **TruncatedIcosidodecahedron**（RS=100% > RZ=66.7%），
これは最大の面数（62面）と複雑な面種混合（四角形・六角形・十角形）を持つ多面体で，
緯度帯が不均一なため max $\varphi$（RS）と max $z$（RZ）の選択が乖離する．

#### R2 逆転現象の消滅

旧版では歪み立方体（Snub Cube）で R2=35.0% > RS=16.7% という「R2 逆転現象」が観測されていた．
改訂版では **R2=0%，RS=RZ=36.7%** となり逆転は消滅した．
原因：局所参照（$\mathbf{c}_{k-1}$）では左条件の対称等変性が崩れ，一部のペアで R2 が誤って左候補を得ていた．
大域参照（$\mathbf{c}_1$）により均一な左条件が適用されると R2 の loxodrome 選択は Snub Cube の三角形密度に対応できず 0% となる．
**R2 逆転現象は 4D の 16胞体のみで観測される現象**であることが明らかになった．

#### 成功率の大幅向上

RS/RZ が 100% を達成する多面体が 0 種（旧）→ 4 種（新）に増加：
- TruncatedOctahedron（旧 75.0% → 新 100%）
- TruncatedCuboctahedron（旧 0% → 新 100%）
- TruncatedIcosahedron（旧 80.6% → 新 100%）
- TruncatedIcosidodecahedron（RS のみ，旧 0% → 新 100%）

これらは面推移的でないため等変性理論による必然ではなく，アルゴリズムの改良による実質的な性能向上である．

#### 0% の多面体

全規則で 0% の多面体は 8種（旧）→ **7種**（新）に減少：
TruncatedCuboctahedron が新たに 100% を達成したため．

---

---

## 13. 作業履歴と再開時の状態（2026-05-07）

### 13.1 このセッションで完了した作業

| 日付 | 作業内容 |
|------|---------|
| 2026-05-06 | アルゴリズム改訂（局所 c_{k-1} 参照 → 大域 c_1 参照，単一候補優先廃止，eps 導入）。`peeling3DLoxo.m` 更新。 |
| 2026-05-06 | `run_platonic_updated.m` を再実行 → `dataPlatonic.mx` 更新。全結果が 0%/100% で等変性理論と一致。 |
| 2026-05-07 | eps < 0（厳密条件，thresh = +10^-10）での RS/RZ 計算を `run_platonic_strict.m` で実行 → `dataPlatonic_strict.mx` 保存。 |
| 2026-05-07 | eps の3ケース（>0/=0/<0）比較を `run_platonic_eps_comparison.m` で実行・確認。 |
| 2026-05-07 | 南極問題（no fallback 失敗の幾何学的原因）を解析し，Section 11.7 として summary.md に追記。 |
| 2026-05-07 | 規則の呼び名を決定（RS = Spiral rule，RZ = Zonal rule）。Section 11.5 に記録。 |
| 2026-05-07 | `visualize_platonic_nets_v2.m` を新規作成（RS/RZ のみ，標準 + 厳密条件の両方を表示）。 |
| 2026-05-07 | `run_archimedean_faceup.m` を改訂版アルゴリズムで再実行 → `archimedean_faceup_results.mx` 更新。Section 12 に結果を追記。 |
| 2026-05-07 | `paper_draft.tex` の Section 2.4 を "Two Selection Rules" に改訂（R2 除去，R1→RS，R3→RZ 統一）。Platonic/Archimedean 両テーブルを更新。 |
| 2026-05-07 | `visualize_archimedean_nets.m` を新規作成（アルキメデス多面体展開図可視化，RS/RZ，`showR2 = False`）。 |
| 2026-05-07 | Kaino2019 著者名を「海野貴也」→「海野啓明」に修正（`summary.md`・`summary.tex`・`summary_en.tex`）。 |

### 13.2 現在の主要データファイルの状態

| ファイル | 内容 | 最終更新 |
|----------|------|---------|
| `dataPlatonic.mx` | 正多面体 R1/R2/R3 × fallback あり/なし，eps > 0 | 2026-05-06 |
| `dataPlatonic_strict.mx` | 正多面体 RS/RZ × fallback あり/なし，eps < 0 | 2026-05-07 |
| `archimedean_faceup_results.mx` | アルキメデス13種 R1/R2/R3，with fallback のみ | 2026-05-07（再計算） |
| `ans4D.mx` | 4D正多胞体 5種の計算結果 | 保存済み（日付不明） |

### 13.3 次回再開時の優先タスク

1. **`peeling4Df3.m` のテスト**（xyz 3成分行列式左条件）：5胞体・8胞体で現行版 `peeling4Df.m` と比較
2. **600胞体の計算完了**：計算量が大きいため分割実行が必要
3. **STL 570枚の valid net 比率集計**：`unfold3DExport.m` の SAT 判定結果を集計
4. **`peeling1back.m` との成功率比較**：1段バックトラック拡張の効果測定

---

## 参考文献

1. T. Yoshino, S. Chaidee, P. Sawatdithep, *Apple-Peel Unfolding of Polyhedra*, `applePeelingV01.pdf`.
2. 海野啓明, *正多胞体の展開図*, 対称性学会 金沢大会 (2019), `Kaino2019_SymmetryCongressKanazawa.pdf`.
