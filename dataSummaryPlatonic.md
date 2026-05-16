# Platonic Solids 計算データの現状

## 保存済み .mx ファイル

| ファイル | 内容 |
|----------|------|
| `dataPlatonic.mx` | 正多面体 5種 × R1/R2/R3 × fallback あり/なし の全ペア計算結果（`run_platonic_updated.m` 実行後に生成） |
| `ans4D.mx` | 4D正多胞体（5・8・16・24・120胞体）の計算結果 |
| `archimedean_results.mx` | アルキメデス13種の旧版結果（Face-up なし・単一候補優先なし） |
| `archimedean_faceup_results.mx` | アルキメデス13種の再計算結果（Face-up あり・ただし単一候補優先なし・R3タイブレーク修正前） |

## dataPlatonic.mx のデータ構造

`run_platonic_updated.m` を実行すると変数 `platonicData` が生成・保存される．

```mathematica
(* 読み込み *)
Get["dataPlatonic.mx"]   (* または絶対パスで *)

(* アクセス例 *)
platonicData["Dodecahedron"]["R3"]["withFallback"]
(* -> {{top, f2, order, success}, ...} のリスト *)
```

### 各エントリの形式

```
{top, f2, order, success}
```

| フィールド | 型 | 説明 |
|------------|-----|------|
| `top` | Integer | 開始面のインデックス |
| `f2` | Integer | 2番目の面のインデックス（実面番号） |
| `order` | List[Integer] | ペーリング順序；成功時は全面数分，失敗時は詰まった面まで |
| `success` | True / False | 全面をペーリングできたか |

### トップレベル構造

```mathematica
platonicData[name] = <|
  "nFaces"     -> Integer,      (* 面数 *)
  "totalPairs" -> Integer,      (* (top, f2) ペアの総数 *)
  "R1" -> <|
    "withFallback" -> { {top, f2, order, success}, ... },
    "noFallback"   -> { {top, f2, order, success}, ... }
  |>,
  "R2" -> <| ... |>,
  "R3" -> <| ... |>
|>
```

`name` は `"Tetrahedron"`, `"Cube"`, `"Octahedron"`, `"Dodecahedron"`, `"Icosahedron"` のいずれか．

### 使用例

```mathematica
(* 成功したものだけ取り出す *)
success = Select[platonicData["Dodecahedron"]["R1"]["withFallback"], #[[4]] &];

(* 失敗したものだけ取り出す（order は途中まで） *)
failed = Select[platonicData["Dodecahedron"]["R2"]["withFallback"], !#[[4]] &];

(* 成功率 *)
all = platonicData["Cube"]["R3"]["noFallback"];
N[Length[Select[all, #[[4]] &]] / Length[all]]
```

## peeling3DLoxo.m の役割

`run_platonic_updated.m` の冒頭で `Get["peeling3DLoxo.m"]` を呼んでいる．
これは3D版アルゴリズムの**コアライブラリ**であり，以下の関数を提供する．

| 関数 | 役割 |
|------|------|
| `p3lExtractPolyhedron[name]` | `PolyhedronData` から頂点・面データを取得 |
| `p3lBuildAdj[faces]` | 面隣接リストを構築 |
| `p3lAlignTopToZ[vers, faces, top]` | 開始面 (top) の重心を +z 軸に揃える回転 |
| `p3lFaceCentroid[fi, vers, faces]` | 面 fi の重心を計算 |
| `p3lDarbouxFrame[cPrev, cCurr]` | Darboux フレーム `{nHat, fHat, lHat}` を構築 |
| `p3lAngleFromForward[fHat, lHat, cj]` | 候補面 cj の前進方向からの角度 φ を計算 |
| `p3lPeelPair[cc, adj, nf, top, f2, rule]` | 1ペア (top, f2) のペーリングを実行（fallback あり，規則 R1/R2/R3） |
| `p3lRunAll[vers, faces, rule]` | 全 (top, f2) ペアを実行 |
| `p3lFromBuiltin[name, rule]` | 多面体名と規則を指定して実行・結果を表示 |
| `p3lCompareRules[]` | 全正多面体で R2・R3 を比較 |

`run_platonic_updated.m` 自体は：
1. `peeling3DLoxo.m` を読み込んで上記関数を使う
2. 追加で **fallback なし版**（`p3lPeelPairNoFB2`・`p3lRunAllNoFB2`）をローカルに定義
3. 5種 × 3規則 × fallback あり/なし の計算を実行し，結果を `dataPlatonic.mx` に保存する

という**実行スクリプト**の位置づけである．

---

## 正多面体の計算スクリプト

| スクリプト | 内容 |
|------------|------|
| `run_platonic_updated.m` | 5種 × R1/R2/R3 × fallback あり/なし の結果を計算し `dataPlatonic.mx` に保存 |
| `run_platonic_fallback.m` | fallback あり/なし の比較を Print するだけ（DumpSave なし・旧版） |

## 計算結果（2026-05-06 版・アルゴリズム改訂後）

アルゴリズム：大域 c_1 参照・`eps = 10^-10`・単一候補優先なし（`run_platonic_updated.m` 実行）

| 多面体 | ペア数 | R1(w) | R1(n) | R2(w) | R2(n) | R3(w) | R3(n) |
|--------|:------:|------:|------:|------:|------:|------:|------:|
| Tetrahedron | 12 | 12 | 12 | 12 | 12 | 12 | 12 |
| Cube | 24 | 24 | 24 | 24 | 0 | 24 | 24 |
| Octahedron | 24 | 24 | 24 | 24 | 0 | 24 | 24 |
| Dodecahedron | 60 | 60 | 60 | 60 | 0 | 60 | 60 |
| Icosahedron | 60 | 60 | 60 | 0 | 0 | 60 | 60 |

w = with fallback，n = no fallback

全結果が 0% か 100%—面推移性による等変性理論と一致。

**備考**：Cube の計算中に `Det::luc` 警告（数値的悪条件行列）が出るが，結果には影響なし。

## 再計算が必要なもの

- 正多面体：再計算済み（結果は `dataPlatonic.mx` に保存，2026-05-06）
- アルキメデス13種：`run_archimedean_faceup.m` を再実行（現行 `peeling3DLoxo.m` は大域 c_1 参照・eps 導入・単一候補優先なし に改訂済みのため，旧結果は無効）

## 可視化（`visualize_platonic_nets.m`）

```mathematica
Get["/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/visualize_platonic_nets.m"]
```

`dataPlatonic.mx` を読み込んで展開図を描画する（内部での再計算はなし）．

### 表示内容

- 5種の正多面体 × R1/R2/R3 × **withFallback / noFallback** を並べて比較表示
- 成功展開図：薄灰色背景・青→緑→赤グラデーション（ペーリング順）
- 失敗展開図：薄灰色背景・× マーク・「配置面数/全面数」ラベル

### 設定パラメータ（スクリプト冒頭で変更可）

| パラメータ | デフォルト | 説明 |
|-----------|:----------:|------|
| `netSize` | 120 | 各展開図の画像サイズ (px) |
| `maxShowOK` | `Infinity` | 成功展開図の表示上限 |
| `maxShowFail` | `Infinity` | 失敗展開図の表示上限 |
| `deduplicate` | `True` | `True`: 同じ peeling order を1件に集約，`False`: 全ペア表示 |

### 実装上の注意

- `p3lUnfoldNet` は面ローカル座標辞書（`fp[k, v]` DownValue）を使用．同じ頂点でも面ごとに異なる2D位置を持てる．
- 外向き法線の保証：`Cross[e1,e2]` の符号を面重心との内積で検証し，内向きなら反転する．
- 重複排除は `DeleteDuplicatesBy[results, #[[3]] &]`（order リストが同一のものを1件に集約）．
- Mathematica の Unicode エスケープは `\:XXXX` 形式を使用（`\uXXXX` は不可）．
