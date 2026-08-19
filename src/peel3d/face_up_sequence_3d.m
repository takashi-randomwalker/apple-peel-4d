(* ================================================================
   face_up_sequence_3d.m
   ----------------------------------------------------------------
   3次元多面体の Apple-Peel 展開について，
   face-up 回転アニメーション制作に必要な情報を書き出す．

   主目的:
     (1) 展開した順序を表す「面番号」のリスト   faceOrder
     (2) 展開に利用した「辺番号」の順序リスト   hingeEdgeIds
   加えて，アニメ制作に有用な補助データ（回転行列・二面角・
   ヒンジ端点・face-up 前後の座標）も同梱する．

   依存: peeling3DLoxo.m（同ディレクトリ）

   使い方:
     Get["face_up_sequence_3d.m"];

     res = p3lPeelSequence["Dodecahedron", 1, "maxz"];
     res["faceOrder"]      (* -> {1, 4, 9, ...}                 *)
     res["hingeEdgeIds"]   (* -> {17, 3, 22, ...}  (長さ nf-1)  *)

     (* すべての多面体・規則で一覧を表示 *)
     p3lPrintSequence["Icosahedron", 1, "maxz"];
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
If[scriptDir === "", scriptDir = NotebookDirectory[]];
Get[scriptDir <> "peeling3DLoxo.m"];


(* ----------------------------------------------------------------
   回転行列 (Rodrigues) — from を to に重ねる proper rotation
   （gen_exampleOfPeeling_v2.m から移植）
   ---------------------------------------------------------------- *)
p3lRotationMatrix3D[from_, to_] :=
  Module[{v, c, s, vx, perp},
    v = Cross[from, to]; c = from . to; s = Norm[v];
    Which[
      s < 10^-10 && c > 0, IdentityMatrix[3],
      s < 10^-10 && c <= 0,
        perp = If[Abs[from[[1]]] < 0.9,
                  Normalize[Cross[from, {1., 0., 0.}]],
                  Normalize[Cross[from, {0., 1., 0.}]]];
        2 Outer[Times, perp, perp] - IdentityMatrix[3],
      True,
        vx = {{0, -v[[3]], v[[2]]}, {v[[3]], 0, -v[[1]]}, {-v[[2]], v[[1]], 0}};
        IdentityMatrix[3] + vx + vx . vx (1 - c) / s^2
    ]];


(* ----------------------------------------------------------------
   稜線リストの構築
     戻り値: {edgeList, edgeId}
       edgeList : {{va,vb}, ...}  各稜線をソート済み頂点ペアで（1-indexed）
       edgeId   : Association  Sort[{va,vb}] -> 稜線番号（1-indexed）
   ---------------------------------------------------------------- *)
p3lEdgeList[faces_] :=
  Module[{edges, edgeList, edgeId},
    edges = Flatten[
      Table[Sort /@ Partition[faces[[i]], 2, 1, {1, 1}], {i, Length[faces]}],
      1];
    edgeList = Union[edges];
    edgeId   = AssociationThread[edgeList -> Range[Length[edgeList]]];
    {edgeList, edgeId}
  ];


(* ----------------------------------------------------------------
   1本の (top, f2) ペアについて face-up + 展開シーケンスを返す
     name : PolyhedronData 名 or {vers,faces}
     top  : 開始面 F1（1-indexed）
     rule : "maxz"(RZ) / "maxphi"(RS) / "loxo"(R2)
     f2   : 2番目の面（Automatic なら top の最初の成功隣接面を採用）
   戻り値: Association（下記キー参照）
   ---------------------------------------------------------------- *)
p3lPeelSequence[name_, top_, rule_, f2_: Automatic] :=
  Module[{vers, faces, nf, adj, versC, cc1, R, lv, cc,
          edgeList, edgeId, f2use, order, complete,
          f2Cands, results, ok, pick,
          hingePairs, hingeIds, hingeVerts, angles,
          verForFaces, normals, i, a, b, sh, n1, n2, cen1},

    (* --- データ取得 --- *)
    {vers, faces} = If[StringQ[name],
      p3lExtractPolyhedron[name],
      name (* {vers,faces} を直接渡した場合 *)
    ];
    nf  = Length[faces];
    adj = p3lBuildAdj[faces];

    (* --- 重心移動（centered）と face-up 回転行列 R --- *)
    versC = N[vers - ConstantArray[Mean[vers], Length[vers]]];
    cc1   = Mean[versC[[ faces[[top]] ]]];
    R     = If[Norm[cc1] < 10^-10,
              IdentityMatrix[3],
              p3lRotationMatrix3D[Normalize[cc1], {0., 0., 1.}]];
    (* アルゴリズムが用いる face-up 座標（等変性のため p3lAlignTopToZ を使用） *)
    lv = p3lAlignTopToZ[vers, faces, top];
    cc = Table[p3lFaceCentroid[i, lv, faces], {i, nf}];

    (* --- 稜線番号 --- *)
    {edgeList, edgeId} = p3lEdgeList[faces];

    (* --- f2 の決定と展開順序の計算 --- *)
    If[f2 === Automatic,
      (* top の隣接面のうち，最初に完全展開できるものを選ぶ *)
      f2Cands = adj[[top]];
      pick = Null;
      Do[
        {order, complete} = p3lPeelPair[cc, adj, nf, top, f2Cands[[fi]], rule];
        If[complete, pick = f2Cands[[fi]]; Break[]],
        {fi, Length[f2Cands]}];
      If[pick === Null,
        (* どれも完全展開しなければ最初の隣接面で（不完全でも）返す *)
        f2use = f2Cands[[1]];
        {order, complete} = p3lPeelPair[cc, adj, nf, top, f2use, rule],
        f2use = pick;
        {order, complete} = p3lPeelPair[cc, adj, nf, top, f2use, rule]
      ],
      f2use = f2;
      {order, complete} = p3lPeelPair[cc, adj, nf, top, f2use, rule]
    ];

    (* --- ヒンジ辺（連続する2面の共有稜線）と二面角 --- *)
    verForFaces = Map[lv[[#]] &, faces, {2}];
    normals = Table[
      With[{verts = verForFaces[[i]]},
        With[{nn = Normalize[Cross[verts[[2]] - verts[[1]],
                                   verts[[3]] - verts[[1]]]]},
          If[nn . Mean[verts] < 0, -nn, nn]]],
      {i, nf}];

    hingePairs = {}; hingeIds = {}; hingeVerts = {}; angles = {};
    Do[
      a  = order[[i - 1]];
      b  = order[[i]];
      sh = Intersection[faces[[a]], faces[[b]]];   (* 共有稜線 = 2頂点 *)
      If[Length[sh] < 2,
        (* 隣接でない（本来起きない）: プレースホルダ *)
        AppendTo[hingePairs, sh]; AppendTo[hingeIds, 0];
        AppendTo[hingeVerts, {}]; AppendTo[angles, 0.];
        Continue[]];
      sh = Sort[sh[[1 ;; 2]]];
      AppendTo[hingePairs, sh];
      AppendTo[hingeIds, edgeId[sh]];
      AppendTo[hingeVerts, lv[[sh]]];             (* face-up 座標での端点 *)
      n1 = normals[[a]]; n2 = normals[[b]];
      AppendTo[angles, ArcCos[Clip[n1 . n2, {-1., 1.}]]],
      {i, 2, Length[order]}
    ];

    <|
      "name"             -> If[StringQ[name], name, "custom"],
      "rule"             -> rule,
      "top"              -> top,
      "f2"               -> f2use,
      "nFaces"           -> nf,
      "complete"         -> complete,
      (* ★ 主要出力 *)
      "faceOrder"        -> order,              (* 面番号の展開順（長さ nf） *)
      "hingeEdgeIds"     -> hingeIds,           (* 辺番号の順序（長さ nf-1） *)
      (* 補助（アニメ用） *)
      "hingeVertexPairs" -> hingePairs,         (* 各ヒンジの頂点番号ペア *)
      "hingeEndpoints"   -> hingeVerts,         (* 各ヒンジ端点の face-up 3D 座標 *)
      "dihedralAngles"   -> angles,             (* 各展開の折り角（rad） *)
      "faceUpRotation"   -> R,                  (* face-up 回転行列（centered→up） *)
      "topCentroidDir"   -> If[Norm[cc1] < 10^-10, {0.,0.,1.}, Normalize[cc1]],
      "verticesCentered" -> versC,              (* 重心移動後（face-up 前） *)
      "verticesFaceUp"   -> lv,                 (* face-up 後（アルゴリズム座標） *)
      "faces"            -> faces,
      "edgeList"         -> edgeList,           (* 稜線番号→頂点ペア *)
      "faceCentroids"    -> cc                  (* face-up 後の面重心 *)
    |>
  ];


(* ----------------------------------------------------------------
   結果を読みやすく表示
   ---------------------------------------------------------------- *)
p3lPrintSequence[name_, top_, rule_, f2_: Automatic] :=
  Module[{res},
    res = p3lPeelSequence[name, top, rule, f2];
    Print["== ", res["name"], "  (rule: ", rule,
          ", top: ", res["top"], ", f2: ", res["f2"], ") =="];
    Print["  faces        : ", res["nFaces"],
          "   complete: ", res["complete"]];
    Print["  faceOrder    : ", res["faceOrder"]];
    Print["  hingeEdgeIds : ", res["hingeEdgeIds"]];
    Print["  hingePairs   : ", res["hingeVertexPairs"]];
    Print["  dihedral(deg): ",
          Round[res["dihedralAngles"]/Degree, 0.1]];
    res
  ];


Print["face_up_sequence_3d.m loaded."];
Print["Usage:"];
Print["  res = p3lPeelSequence[\"Dodecahedron\", 1, \"maxz\"]  -- 情報を Association で取得"];
Print["  p3lPrintSequence[\"Icosahedron\", 1, \"maxz\"]        -- 見やすく表示"];
Print["  keys: faceOrder, hingeEdgeIds, hingeVertexPairs, dihedralAngles,"];
Print["        faceUpRotation, verticesCentered, verticesFaceUp, edgeList, ..."];
