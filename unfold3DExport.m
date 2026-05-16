(* ================================================================
   unfold3DExport.m
   ================================================================
   Apple-peel peeling の結果 (results4D) を読み込み，
   連結順に従って各胞を独立した正多面体として3次元空間に配置し，
   STLファイルとして出力する．

   各胞は独立した頂点座標を持つ．連続する胞間では共有面の
   頂点座標が完全に一致する．

   前提:
     - peeling4Df.m が読み込まれていること
     - results4D が定義済みであること (4DPeelingAll.nb 実行後)
     - dataDir が _4DData/ フォルダを指していること

   使い方:
     Get[NotebookDirectory[] <> "unfold3DExport.m"];
     processAllResults[];     (* 全多胞体を処理 *)

     (* 個別処理も可能 *)
     processOnePolytope[5];   (* 5胞体のみ *)
   ================================================================ *)

(* === 補助関数 === *)

(* 胞 ci の頂点番号リスト (1-indexed) *)
unfoldCellVerts[ci_, cells_, faces_] :=
  Union @@ faces[[ cells[[ci]] ]];

(* 胞 ci の4D頂点をSVD射影して局所3D座標を得る．
   各胞は3次元超平面上にあるため，この射影は等長写像 *)
unfoldCellLocal3D[ci_, vers_, cells_, faces_] :=
  Module[{vi, P, c, Pc, svd, basis, local3d},
    vi = unfoldCellVerts[ci, cells, faces];
    P = N[vers[[vi]]];
    c = Mean[P];
    Pc = P - ConstantArray[c, Length[P]];
    svd = SingularValueDecomposition[Pc];
    (* svd[[3]] = V 行列, 最初の3列が3D超平面の基底 *)
    basis = svd[[3]][[All, 1 ;; 3]];
    local3d = Pc . basis;
    AssociationThread[vi -> Table[local3d[[i]], {i, Length[vi]}]]
  ];

(* 3点以上から面の法線ベクトルを計算 *)
unfoldFaceNormal[pts_] :=
  Module[{n},
    n = Cross[pts[[2]] - pts[[1]], pts[[3]] - pts[[1]]];
    If[Norm[n] > 10^-10, Normalize[N[n]], {0., 0., 1.}]
  ];

(* Procrustes 整列 (スケールなし): R.src + t -> tgt *)
unfoldProcrustes[srcPts_, tgtPts_] :=
  Module[{cS, cT, Sc, Tc, H, svd, u, v, d, R, t},
    cS = Mean[srcPts]; cT = Mean[tgtPts];
    Sc = srcPts - ConstantArray[cS, Length[srcPts]];
    Tc = tgtPts - ConstantArray[cT, Length[tgtPts]];
    H = Transpose[Sc] . Tc;
    svd = SingularValueDecomposition[H];
    u = svd[[1]]; v = svd[[3]];
    d = Det[v . Transpose[u]];
    R = v . DiagonalMatrix[{1., 1., d}] . Transpose[u];
    t = cT - R . cS;
    {R, t}
  ];

(* === メイン展開関数 === *)

(* 連結順 order に従い，各胞を独立した正多面体として3D配置する．
   戻り値: {emb1, emb2, ...}  各 emb は Association[頂点番号 -> {x,y,z}]
   cellLocal3DsOpt: オプション引数．Length[cells] の事前計算済み局所座標リスト．
   渡すと unfoldCellLocal3D の重複 SVD を省略できる（大型胞体で高速化）． *)
unfoldTo3DCore[cells_, faces_, order_, getLocal3D_] :=
  Module[{nc, cellEmbs, prevCI, currCI, prevEmb,
          sharedFI, sv, localCurr, currVI,
          A, B, R, t, currEmb,
          faceCtr, nFace, prevCenter, prevSign,
          newCenter, newSign, fcN},

    nc = Length[order];
    (* 事前確保: AppendTo の O(n^2) コピーを回避 *)
    cellEmbs = Table[Null, {nc}];
    cellEmbs[[1]] = getLocal3D[order[[1]]];

    Do[
      prevCI  = order[[k]];
      currCI  = order[[k + 1]];
      prevEmb = cellEmbs[[k]];

      (* 共通する面番号 *)
      sharedFI = Sort[Intersection[cells[[prevCI]], cells[[currCI]]]];

      If[Length[sharedFI] == 0,
        (* 共有面なし（通常は起こらない）: 局所配置のみ *)
        cellEmbs[[k + 1]] = getLocal3D[currCI],

        (* --- 共有面あり: Procrustes整列 + 必要なら鏡映 --- *)

        (* 共有面の頂点集合 *)
        sv = Sort[Union @@ faces[[sharedFI]]];

        (* 新しい胞の局所3D座標 *)
        localCurr = getLocal3D[currCI];
        currVI = unfoldCellVerts[currCI, cells, faces];

        (* Procrustes: localCurr の共有面頂点 -> prevEmb の共有面頂点 *)
        A = localCurr /@ sv;
        B = prevEmb /@ sv;
        {R, t} = unfoldProcrustes[A, B];

        (* 全頂点を変換 *)
        currEmb = AssociationThread[
          currVI -> Table[R . localCurr[v] + t, {v, currVI}]
        ];

        (* 面の法線と前の胞の中心の符号を調べる *)
        faceCtr = Mean[B];
        nFace = unfoldFaceNormal[B];

        prevCenter = Mean[prevEmb /@ unfoldCellVerts[prevCI, cells, faces]];
        prevSign = (prevCenter - faceCtr) . nFace;

        newCenter = Mean[currEmb /@ currVI];
        newSign = (newCenter - faceCtr) . nFace;

        (* 同じ側にあれば面を通して鏡映 *)
        If[prevSign * newSign > 0,
          fcN = faceCtr . nFace;
          currEmb = AssociationThread[
            currVI -> Table[
              currEmb[v] - 2 (currEmb[v] . nFace - fcN) nFace,
              {v, currVI}
            ]
          ];
        ];

        cellEmbs[[k + 1]] = currEmb
      ],
      {k, 1, nc - 1}
    ];

    cellEmbs
  ];

(* 標準版: vers から毎回 unfoldCellLocal3D を呼ぶ（後方互換） *)
unfoldTo3D[vers_, faces_, cells_, order_] :=
  unfoldTo3DCore[cells, faces, order,
    Function[ci, unfoldCellLocal3D[ci, vers, cells, faces]]];

(* 高速版: 事前計算済み cellLocal3Ds（Length[cells] のリスト）を受け取る
   大型胞体（120胞体など）で 1 order あたりの SVD を省略できる *)
unfoldTo3DFast[cells_, faces_, order_, cellLocal3Ds_] :=
  unfoldTo3DCore[cells, faces, order,
    Function[ci, cellLocal3Ds[[ci]]]];

(* === 3D重なり判定 (Separating Axis Theorem) === *)

(* 胞のSAT判定用データを事前計算 *)
unfoldCellSATData[ci_, cells_, faces_, emb_] :=
  Module[{vi, pts, normals, edgeDirs, edgeSet, fi, vlist, p,
          nm, nn, f, a, b, d, dn},
    vi = unfoldCellVerts[ci, cells, faces];
    pts = N[emb /@ vi];
    (* 面法線 *)
    normals = {};
    Do[
      vlist = faces[[fi]];
      p = emb /@ vlist;
      If[Length[p] >= 3,
        nm = Cross[p[[2]] - p[[1]], p[[3]] - p[[1]]];
        nn = Norm[nm];
        If[nn > 10^-10, AppendTo[normals, nm / nn]]
      ],
      {fi, cells[[ci]]}
    ];
    (* 辺方向ベクトル (正規化) *)
    edgeSet = {};
    Do[
      f = faces[[fi]];
      Do[
        b = Mod[a, Length[f]] + 1;
        AppendTo[edgeSet, Sort[{f[[a]], f[[b]]}]],
        {a, Length[f]}
      ],
      {fi, cells[[ci]]}
    ];
    edgeSet = DeleteDuplicates[edgeSet];
    edgeDirs = {};
    Do[
      d = emb[edgeSet[[i, 2]]] - emb[edgeSet[[i, 1]]];
      dn = Norm[d];
      If[dn > 10^-10, AppendTo[edgeDirs, d / dn]],
      {i, Length[edgeSet]}
    ];
    Association[
      "pts" -> pts,
      "normals" -> normals,
      "edges" -> edgeDirs,
      "bboxMin" -> Min /@ Transpose[pts],
      "bboxMax" -> Max /@ Transpose[pts]
    ]
  ];

(* SAT で2つの凸多面体の重なりを判定
   eps > 0 により共有面での接触は重なりと見なさない
   最適化: 面法線を先に検査し，分離確定なら辺×辺を生成しない（Throw/Catch）*)
satOverlap[cdA_, cdB_, eps_: 10^-6] :=
  Module[{ptsA, ptsB, projA, projB, c, cn},
    ptsA = cdA["pts"];
    ptsB = cdB["pts"];
    Catch[
      (* Phase 1: 面法線のみ（計 |nA|+|nB| 本，典型的に 24 本） *)
      Do[
        projA = ptsA . axis;
        projB = ptsB . axis;
        If[Max[projA] <= Min[projB] + eps || Max[projB] <= Min[projA] + eps,
          Throw[False]  (* False = 分離済み = 重なりなし *)
        ],
        {axis, Join[cdA["normals"], cdB["normals"]]}
      ];
      (* Phase 2: 辺×辺をオンデマンドで生成・検査（面法線で分離できなかった場合のみ） *)
      Do[
        Do[
          c = Cross[ea, eb]; cn = Norm[c];
          If[cn > 10^-10,
            projA = ptsA . (c/cn);
            projB = ptsB . (c/cn);
            If[Max[projA] <= Min[projB] + eps || Max[projB] <= Min[projA] + eps,
              Throw[False]
            ]
          ],
          {eb, cdB["edges"]}
        ],
        {ea, cdA["edges"]}
      ];
      True  (* 全軸で分離なし → 重なり *)
    ]
  ];

(* 全胞ペアの重なりを判定 *)
checkOverlaps[order_, cells_, faces_, cellEmbs_, eps_: 10^-6] :=
  Module[{n, cds, overlaps, i, j},
    n = Length[order];
    (* 事前計算 *)
    cds = Table[
      unfoldCellSATData[order[[idx]], cells, faces, cellEmbs[[idx]]],
      {idx, n}
    ];
    (* 全ペア検査 (AABB事前フィルタ付き) *)
    overlaps = {};
    Do[
      Do[
        (* AABB事前フィルタ *)
        If[Or @@ Thread[cds[[i]]["bboxMax"] < cds[[j]]["bboxMin"] - eps],
          Continue[]];
        If[Or @@ Thread[cds[[j]]["bboxMax"] < cds[[i]]["bboxMin"] - eps],
          Continue[]];
        (* SAT詳細判定 *)
        If[satOverlap[cds[[i]], cds[[j]], eps],
          AppendTo[overlaps, {i, j}]
        ],
        {j, i + 1, n}
      ],
      {i, 1, n - 1}
    ];
    overlaps
  ];

(* === 三角形の収集 === *)

unfoldTriangles[order_, cells_, faces_, cellEmbs_] :=
  Module[{tris = {}, ci, emb, vlist, pts},
    Do[
      ci = order[[idx]];
      emb = cellEmbs[[idx]];
      Do[
        vlist = faces[[fi]];
        pts = emb /@ vlist;
        (* 扇形三角形分割 *)
        Do[
          AppendTo[tris, {pts[[1]], pts[[j]], pts[[j + 1]]}],
          {j, 2, Length[pts] - 1}
        ],
        {fi, cells[[ci]]}
      ],
      {idx, Length[order]}
    ];
    tris
  ];

(* === 可視化 === *)

showUnfolding[vers_, faces_, cells_, order_] :=
  Module[{cellEmbs, nCells, colors, polys, ci, emb, col},
    cellEmbs = unfoldTo3D[vers, faces, cells, order];
    nCells = Length[order];
    colors = Table[ColorData[97][i], {i, nCells}];
    polys = {};
    Do[
      ci = order[[idx]];
      emb = cellEmbs[[idx]];
      col = colors[[idx]];
      Do[
        AppendTo[polys,
          {EdgeForm[{Thin, Darker[col, 0.3]}],
           FaceForm[Directive[col, Opacity[0.4]]],
           Polygon[emb /@ faces[[fi]]]}
        ],
        {fi, cells[[ci]]}
      ],
      {idx, nCells}
    ];
    Graphics3D[polys, Boxed -> False, Lighting -> "Neutral"]
  ];

(* === STL出力 === *)

(* 1つの連結順を STL に出力 *)
exportUnfoldSTL[vers_, faces_, cells_, order_, filename_] :=
  Module[{cellEmbs, tris, gr},
    cellEmbs = unfoldTo3D[vers, faces, cells, order];
    tris = unfoldTriangles[order, cells, faces, cellEmbs];
    gr = Graphics3D[{EdgeForm[], Polygon /@ tris}];
    Export[filename, gr, "STL"];
    Print["  Exported: ", filename,
      "  (", Length[order], " cells, ", Length[tris], " triangles)"];
    filename
  ];

(* === 全結果の一括処理 === *)

processOnePolytope[n_] := processOnePolytope[n,
  If[ValueQ[dataDir], dataDir, NotebookDirectory[] <> "_4DData/"],
  If[ValueQ[outDir],  outDir,  NotebookDirectory[]]
];

processOnePolytope[n_, dataDir0_, outDir0_] :=
  Module[{raw, vers, faces, cells, allOrders, ord, success,
          count, fname, cellEmbs, overlaps, nOK, nOV,
          overlapResults},

    Print["=== ", n, "-cell ==="];

    (* 元データ読み込み: {vers, edgs, neis, faces, fneis, cells} *)
    raw = Get[dataDir0 <> "f" <> ToString[n] <> ".m"];
    vers  = raw[[1]];
    faces = raw[[4]];
    cells = raw[[6]];

    (* results4D から成功した全連結順を収集 *)
    allOrders = {};
    Do[
      Do[
        ord     = results4D[n]["results"][[top]][[f2]][[2]];
        success = results4D[n]["results"][[top]][[f2]][[3]];
        If[success, AppendTo[allOrders, ord]],
        {f2, Length[results4D[n]["results"][[top]]]}
      ],
      {top, Length[results4D[n]["results"]]}
    ];

    (* 重複除去 *)
    allOrders = DeleteDuplicates[allOrders];
    Print["  Unique successful orders: ", Length[allOrders]];

    (* 各連結順で展開 → 重なり判定 → STL出力 *)
    count = 0; nOK = 0; nOV = 0;
    overlapResults = {};

    Do[
      count++;
      cellEmbs = unfoldTo3D[N[vers], faces, cells, allOrders[[i]]];
      overlaps = checkOverlaps[allOrders[[i]], cells, faces, cellEmbs];

      If[Length[overlaps] == 0,
        nOK++;
        Print["  [", count, "/", Length[allOrders], "] OK"],
        nOV++;
        Print["  [", count, "/", Length[allOrders], "] OVERLAP  (",
              Length[overlaps], " pairs)"]
      ];
      AppendTo[overlapResults,
        Association["order" -> allOrders[[i]],
                    "overlaps" -> overlaps,
                    "hasOverlap" -> (Length[overlaps] > 0)]
      ];

      fname = outDir0 <> "unfold_" <> ToString[n] <> "cell_" <>
              ToString[count] <> ".stl";
      exportUnfoldSTL[N[vers], faces, cells, allOrders[[i]], fname],
      {i, Length[allOrders]}
    ];

    Print["  Summary: ", nOK, " OK, ", nOV, " OVERLAP"];
    Association["n" -> n, "nOK" -> nOK, "nOverlap" -> nOV,
                "total" -> Length[allOrders],
                "details" -> overlapResults]
  ];

processAllResults[] := processAllResults[
  If[ValueQ[dataDir], dataDir, NotebookDirectory[] <> "_4DData/"],
  If[ValueQ[outDir],  outDir,  NotebookDirectory[]]
];

processAllResults[dataDir0_, outDir0_] :=
  Module[{nList},
    nList = Keys[results4D];
    Print["Processing polytopes: ", nList];
    Do[
      processOnePolytope[n, dataDir0, outDir0],
      {n, nList}
    ];
    Print["=== All done. ==="];
  ];

Print["unfold3DExport.m loaded."];
Print["Usage:"];
Print["  processAllResults[]          -- export all polytopes"];
Print["  processOnePolytope[8]        -- export 8-cell only"];
Print["  showUnfolding[vers,faces,cells,order]  -- visualize"];
