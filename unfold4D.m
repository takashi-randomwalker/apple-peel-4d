(* ================================================================
   unfold4D.m
   ================================================================
   4D正多胞体の3D展開図可視化
   peeling4Df4.m (グローバル c1-c2 参照 + Rule 3) の結果を用いて
   各セルを共有面を軸に回転させ 3D ネットとして描画する．

   アルゴリズム (unfold3DExport.m の unfoldTo3D による):
     1. C1: SVD 射影で 4D → 3D ローカル座標
     2. C_k: SVD 射影 → Procrustes 整列 (共有面の頂点を C_{k-1} に合わせる)
     3. 新セルが前セルと同じ側にあれば共有面について鏡映して「外側」に展開

   使い方 (Mathematica Notebook):
     SetDirectory["/path/to/260324Peeling4D"];
     Get["unfold4D.m"]

   設定パラメータは下の「設定」節で変更できる．
   ================================================================ *)

(* ---- ディレクトリ設定 ---- *)
scriptDir = If[ValueQ[$InputFileName], DirectoryName[$InputFileName],
               NotebookDirectory[]];
dataDir = scriptDir <> "_4DData/";

(* ---- 依存ファイル読み込み ---- *)
Get[scriptDir <> "unfold3DExport.m"];  (* unfoldTo3D など *)
If[FileExistsQ[scriptDir <> "ans4DGlobal_v4.mx"],
  Get[scriptDir <> "ans4DGlobal_v4.mx"],   (* 全6胞体 v4 版（k>=3 Det=0 fix 適用済み） *)
  Get[scriptDir <> "ans4DGlobal.mx"]        (* フォールバック：v4 未生成時 *)
];

(* ================================================================
   設定
   ================================================================ *)
If[!ValueQ[netImageSize],   netImageSize   = 500];   (* Graphics3D の画像サイズ (px) *)
If[!ValueQ[faceOpacity],   faceOpacity    = 0.80];  (* 面の不透明度 *)
If[!ValueQ[viewPt],        viewPt         = {2.5, -2.5, 2.0}]; (* 視点 *)
If[!ValueQ[deduplicateGeo],deduplicateGeo = False];  (* True: 3D形状で重複除去（対称群で等価な net は1件） *)

(* ================================================================
   補助関数
   ================================================================ *)

(* --- 多胞体データ読み込み --- *)
load4DData[n_] :=
  Module[{raw},
    raw = Get[dataDir <> "f" <> ToString[n] <> ".m"];
    Association["vers" -> raw[[1]], "faces" -> raw[[4]], "cells" -> raw[[6]]]
  ];

(* --- Apple-Peel グラデーション色: 濃灰(C1, 0.5) → 淡灰(末尾, 0.9) --- *)
apColor[k_, nTotal_] :=
        GrayLevel[0.5+ 0.4 * (k - 1) / Max[nTotal - 1, 1]];

(* --- 3D 展開図描画 --- *)
(* vers:  face-up 回転後の 4D 頂点座標 (peeling4Df4 の localVers)  *)
(* faces: {頂点インデックスリスト, ...}                              *)
(* cells: {面インデックスリスト, ...}                                *)
(* order: ペーリング順序 (セルインデックスのリスト)                  *)
showNet4D[vers_, faces_, cells_, order_, sz_: 500] :=
  Module[{cellEmbs, nCells, prims, ci, emb, col},
    cellEmbs = unfoldTo3D[vers, faces, cells, order];
    nCells   = Length[order];
    prims = Flatten[Table[
      ci  = order[[idx]];
      emb = cellEmbs[[idx]];
      col = apColor[idx, nCells];
      Map[
        Function[fi,
          {EdgeForm[{Thin, Black}],
           FaceForm[Directive[col, Opacity[faceOpacity]]],
           Polygon[emb /@ faces[[fi]]]}
        ],
        cells[[ci]]
      ]
    , {idx, nCells}], 1];
    Graphics3D[
      prims,
      Boxed        -> False,
      ImageSize    -> sz,
      Lighting     -> "Neutral",
      ViewPoint    -> viewPt,
      SphericalRegion -> True
    ]
  ];

(* --- 幾何学的重複除去キー (4D → 3D 展開図)
   セル重心の正規化済みペアワイズ距離行列（ソート済み）をキーとして返す．
   回転・並進・スケール不変なため，対称群で等価な異なる order が同じキーになる．
   cellEmbs: unfoldTo3D の戻り値（セルごとの埋め込み関数のリスト） *)
netGeoKey4D[cellEmbs_, faces_, cells_, order_] :=
  Module[{centroids, c1, shifted, maxD},
    centroids = Table[
      Mean[cellEmbs[[idx]] /@ (Union @@ faces[[ cells[[ order[[idx]] ]] ]])],
      {idx, Length[order]}
    ];
    c1      = centroids[[1]];
    shifted = # - c1 & /@ centroids;
    maxD    = Max[Norm /@ shifted];
    If[maxD < 10^-8, Return[{}]];
    shifted = shifted / maxD;
    Round[Sort@Flatten@Table[
      Norm[shifted[[i]] - shifted[[j]]],
      {i, Length[shifted]}, {j, i + 1, Length[shifted]}
    ], 10^-3]
  ];

(* --- 既計算の cellEmbs から 3D 展開図を描画（unfoldTo3D 再計算不要） --- *)
showNet4DFromEmb[cellEmbs_, faces_, cells_, order_, sz_: 500] :=
  Module[{nCells, prims},
    nCells = Length[order];
    prims  = Flatten[Table[
      Module[{ci = order[[idx]], emb = cellEmbs[[idx]], col = apColor[idx, nCells]},
        Map[Function[fi,
              {EdgeForm[{Thin, Black}],
               FaceForm[Directive[col, Opacity[faceOpacity]]],
               Polygon[emb /@ faces[[fi]]]}],
            cells[[ci]]]
      ],
      {idx, nCells}], 1];
    Graphics3D[prims, Boxed -> False, ImageSize -> sz,
      Lighting -> "Neutral",
      ViewPoint -> viewPt, SphericalRegion -> True]
  ];

(* --- 固定 C1 の全成功ペアを取得 (order で重複除去) --- *)
getPairsForTop[nPoly_, top_] :=
  Module[{res},
    res = Select[results4DGlobal[nPoly]["results"][[top]], Last[#] &];
    DeleteDuplicatesBy[res, #[[2]] &]
  ];

(* --- face-up 回転後の頂点を再計算 (カスタム頂点座標を使う場合) --- *)
recomputeLocalVers[vers_, faces_, cells_, top_] :=
  Module[{wHat, localVers, cc1, angle4D, cellVerts},
    wHat = {0., 0., 0., 1.};
    cellVerts[i_] := Union @@ faces[[ cells[[i]] ]];
    localVers = (# - Mean[vers]) & /@ vers;
    cc1 = Mean[ localVers[[ cellVerts[top] ]] ];
    Which[
      Abs[Normalize[cc1] . wHat - 1] < 0.0001, localVers,
      Abs[Normalize[cc1] . wHat + 1] < 0.0001, -localVers,
      True,
        angle4D = ArcCos[Clip[cc1 . wHat / Norm[cc1], {-1, 1}]];
        RotationMatrix[angle4D, {Normalize[cc1], wHat}] . # & /@ localVers
    ]
  ];

(* --- 多胞体1種分を描画するヘルパー ---                           *)
(* customVers: 頂点座標を上書きする場合に指定 (省略時は保存済み localVers を使用) *)
showPolytopeNets[nPoly_, poly_, top_, label_, col_, customVers_: None] :=
  Module[{pairs, nOrd, n, lv, ord, seenKeys, geoUniq, emb, key},
    Print[Style["━━━ " <> label <> " ━━━", Bold, 16, col]];
    pairs = getPairsForTop[nPoly, top];
    nOrd  = Length[pairs];

    If[deduplicateGeo,

      (* 幾何学的重複除去: unfoldTo3D を事前計算してキーで絞る *)
      seenKeys = <||>;
      geoUniq  = {};
      Do[
        ord = pairs[[i, 2]];
        lv  = If[customVers === None,
                 pairs[[i, 1]],
                 recomputeLocalVers[customVers, poly["faces"], poly["cells"], ord[[1]]]];
        emb = unfoldTo3D[lv, poly["faces"], poly["cells"], ord];
        key = netGeoKey4D[emb, poly["faces"], poly["cells"], ord];
        If[!KeyExistsQ[seenKeys, key],
          seenKeys[key] = True;
          AppendTo[geoUniq, {lv, ord, emb}]
        ],
        {i, nOrd}
      ];
      n = Length[geoUniq];
      Print["  C1 = ", top, "  →  order-unique: ", nOrd,
            ",  geo-unique: ", n];
      Do[
        {lv, ord, emb} = geoUniq[[i]];
        Print[Style["  [#" <> ToString[i] <> "]  order: " <> ToString[ord],
                    11, Darker[col, 0.2]]];
        Print[showNet4DFromEmb[emb, poly["faces"], poly["cells"], ord, netImageSize]],
        {i, n}],

      (* デフォルト: order ベース重複除去のみ *)
      n = nOrd;
      Print["  C1 = ", top, "  →  ", n, " unique successful net(s)."];
      Do[
        ord = pairs[[i, 2]];
        lv  = If[customVers === None,
                 pairs[[i, 1]],
                 recomputeLocalVers[customVers, poly["faces"], poly["cells"], ord[[1]]]];
        Print[Style["  [#" <> ToString[i] <> "]  order: " <> ToString[ord],
                    11, Darker[col, 0.2]]];
        Print[showNet4D[lv, poly["faces"], poly["cells"], ord, netImageSize]],
        {i, n}]
    ]
  ];

(* ================================================================
   固定 C1 インデックス
   ================================================================ *)
fixedTop = 3;

(* ================================================================
   5胞体 (5 cells, each a tetrahedron)
   ================================================================ *)
poly5 = load4DData[5];
showPolytopeNets[5, poly5, fixedTop,
  "5-cell (pentachoron)", RGBColor[0.1, 0.3, 0.7]];

(* ================================================================
   8胞体 (8 cells, each a cube)
   ================================================================ *)
poly8 = load4DData[8];
showPolytopeNets[8, poly8, fixedTop,
  "8-cell (tesseract)", RGBColor[0.1, 0.5, 0.2]];

(* ================================================================
   16胞体 (16 cells, each a tetrahedron)
   ================================================================ *)
poly16 = load4DData[16];
showPolytopeNets[16, poly16, fixedTop,
  "16-cell (hexadecachoron)", RGBColor[0.5, 0.1, 0.5]];

(* ================================================================
   24胞体 (24 cells, each a regular octahedron)
   ================================================================ *)
poly24 = load4DData[24];
showPolytopeNets[24, poly24, fixedTop,
  "24-cell (icositetrachoron)", RGBColor[0.6, 0.35, 0.0]];

(* ================================================================
   120胞体 (120 cells, each a regular dodecahedron)
   ================================================================ *)
poly120 = load4DData[120];
showPolytopeNets[120, poly120, fixedTop,
  "120-cell (hecatonicosachoron)", RGBColor[0.0, 0.4, 0.6]];

Print[Style["unfold4D.m: done.", Bold, GrayLevel[0.4]]];
