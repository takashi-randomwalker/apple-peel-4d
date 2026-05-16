(* ================================================================
   peeling3DLoxo.m
   ================================================================
   3次元多面体における Apple-Peel ペーリング
     規則2 (loxodrome / 最小角度 "loxo")
     規則3 (max z               "maxz")
   の実装と比較．

   右半空間条件（外側から見て CW = 右利きが左方向に剥く）:
     k >= 1 全て: Det[{c_1, c_k, c_j}] <= eps  (eps = 10^-10)

   選択規則:
     RS "maxphi": 右候補から min Det（最も負 = 最大右旋回 = 外側から CW）
     R2 "loxo": 右候補から max Det（最も 0 に近い = loxodrome 近似）
     RZ "maxz": 右候補から max z（緯度優先），タイは min Det
     フォールバック (右候補なし): RS → min Det，R2 → max |Det|，RZ → min z → min Det

   使い方:
     Get["peeling3DLoxo.m"];

     (* Platonic 立体で両規則を比較 *)
     cmpResult = p3lCompareRules[];

     (* 個別実行 *)
     res = p3lFromBuiltin["Dodecahedron", "loxo"];
   ================================================================ *)


(* ================================================================
   Section 1: 補助関数
   ================================================================ *)

(* 面の重心 *)
p3lFaceCentroid[fi_, vers_, faces_] :=
  Mean[N[vers[[ faces[[fi]] ]]]];

(* 面隣接リスト: 辺共有 = 2頂点以上共有 *)
p3lBuildAdj[faces_] :=
  Module[{nf = Length[faces]},
    Table[
      Select[Range[nf],
        Function[j, j != i &&
          Length[Intersection[faces[[i]], faces[[j]]]] >= 2]],
      {i, nf}
    ]
  ];

(* top 面の重心を +z 軸に整列させた頂点座標を返す *)
p3lAlignTopToZ[vers0_, faces_, top_] :=
  Module[{vers, cc1, zHat, u, dot, vVec, a, R3},
    vers = N[vers0 - ConstantArray[Mean[vers0], Length[vers0]]];
    cc1  = Mean[vers[[ faces[[top]] ]]];
    zHat = {0., 0., 1.};
    If[Norm[cc1] < 10^-10, Return[vers]];
    u   = Normalize[cc1];
    dot = Clip[u . zHat, {-1., 1.}];
    Which[
      Abs[dot - 1] <= 10^-4, vers,
      (* Antipodal case: -vers (inversion) has det = -1 and reverses the
         left-half-space condition.  Use a proper 180° rotation (x-axis)
         instead: (x,y,z) -> (x,-y,-z), det = +1. *)
      Abs[dot + 1] <= 10^-4, Map[(# * {1., -1., -1.}) &, vers],
      True,
        vVec = Normalize[zHat - (zHat . u) u];
        a    = ArcCos[dot];
        R3   = IdentityMatrix[3]
               + Sin[a]  (Outer[Times, vVec, u] - Outer[Times, u, vVec])
               + (Cos[a] - 1)(Outer[Times, u, u] + Outer[Times, vVec, vVec]);
        vers . Transpose[R3]
    ]
  ];

(* Darboux フレームを構築
   戻り値: {nHat, fHat, lHat}
     nHat: 球面法線 = Normalize[c_k]
     fHat: 前進方向 (接平面内, c_{k-1}→c_k の投影)
     lHat: 左方向 = nHat × fHat *)
p3lDarbouxFrame[cPrev_, cCurr_] :=
  Module[{nHat, dStep, fHatRaw, fHat, lHat},
    nHat    = If[Norm[cCurr] > 10^-10, Normalize[cCurr], {0.,0.,1.}];
    dStep   = cCurr - cPrev;
    fHatRaw = dStep - (dStep . nHat) nHat;  (* 接平面へ射影 *)
    fHat    = If[Norm[fHatRaw] > 10^-10,
                 Normalize[fHatRaw],
                 (* フォールバック: nHat に垂直な任意ベクトル *)
                 Normalize[{-nHat[[2]], nHat[[1]], 0.} -
                           ({-nHat[[2]], nHat[[1]], 0.} . nHat) nHat]
               ];
    lHat = Cross[nHat, fHat];  (* = nHat × fHat *)
    {nHat, fHat, lHat}
  ];

(* 候補面 j の前進方向 fHat からの角度 φ_j (ラジアン, 正 = 左)
   範囲: (-π, π] *)
(*  nHat = Cross[fHat, lHat] (右手系フレームより nHat = f × L),
    よって Cross[lHat, fHat] = -nHat．
    cj の接平面射影: cj - (cj · nHat) nHat = cj - (cj·(-Cross[lHat,fHat]))(- Cross[lHat,fHat])
    = cj - (cj · nHat) nHat．
    ArcTan[x, y] は Mathematica の atan2(y, x) → (-π, π] を返す． *)
p3lAngleFromForward[fHat_, lHat_, cj_] :=
  Module[{nHat, proj, pn},
    nHat = Cross[fHat, lHat];   (* = normalize(c_k) *)
    proj = cj - (cj . nHat) nHat;  (* 接平面への射影 *)
    pn   = Norm[proj];
    If[pn < 10^-10, Return[0.]];
    proj = proj / pn;
    ArcTan[Clip[proj . fHat, {-1., 1.}],
           Clip[proj . lHat, {-1., 1.}]]
  ];


(* ================================================================
   Section 2: 1ペア (top, f2) のペーリング
   ================================================================ *)

(*  cc    : 面重心リスト (nf × 3)
    adj   : 隣接リスト
    nf    : 面数
    top   : 開始面 (1-indexed)
    f2    : 2番目の面
    rule  : "loxo" または "maxz"
    戻り値: {order, completeBool} *)
p3lPeelPair[cc_, adj_, nf_, top_, f2_, rule_] :=
  Module[{order, inOrder, nbrs, prev, last,
          pending, leftCands, nHat, fHat, lHat,
          dets, absDets, angles, zvals, maxz, minz,
          tied, tDets, next,
          eps},
    eps = 10^-10;

    order   = {top, f2};
    inOrder = <|top -> True, f2 -> True|>;
    nbrs    = Table[
      DeleteCases[DeleteCases[adj[[i]], top], f2],
      {i, nf}
    ];
    prev = top;
    last = f2;

    While[Length[order] < nf,

      pending = Select[nbrs[[last]], !KeyExistsQ[inOrder, #]&];
      nbrs[[last]] = pending;
      If[Length[pending] == 0, Break[]];

      (* 右半空間: Det[{c_1, c_k, c_j}] <= eps
         Cross[c_1, c_k] との内積が負 = 外側から見て CW = 右利きが左方向に剥く方向 *)
      leftCands = Select[pending, Function[j,
        Det[{cc[[top]], cc[[last]], cc[[j]]}] <= eps
      ]];

      next = Which[

        rule === "maxphi" && Length[leftCands] > 0,
          (* RS: min Det（最も負 = 最大右旋回 = 外側から CW） *)
          dets = Det[{cc[[top]], cc[[last]], cc[[#]]}] & /@ leftCands;
          leftCands[[ First[Ordering[dets, 1]] ]],

        rule === "loxo" && Length[leftCands] > 0,
          (* R2: max Det（最も 0 に近い = 最小右旋回 = loxodrome 近似） *)
          dets = Det[{cc[[top]], cc[[last]], cc[[#]]}] & /@ leftCands;
          leftCands[[ First[Ordering[-dets, 1]] ]],

        rule === "maxz" && Length[leftCands] > 0,
          (* RZ: max z → min Det タイブレーク *)
          zvals = cc[[leftCands, 3]];
          maxz  = Max[zvals];
          tied  = Select[leftCands, cc[[#, 3]] >= maxz - 10^-10 &];
          If[Length[tied] == 1, tied[[1]],
            tDets = Det[{cc[[top]], cc[[last]], cc[[#]]}] & /@ tied;
            tied[[ First[Ordering[tDets, 1]] ]]
          ],

        (* フォールバック (右候補なし) *)
        rule === "maxphi" && Length[pending] > 0,
          (* RS フォールバック: Darboux frame で min φ（最大右旋回） *)
          {nHat, fHat, lHat} = p3lDarbouxFrame[cc[[prev]], cc[[last]]];
          angles = p3lAngleFromForward[fHat, lHat, cc[[#]]] & /@ pending;
          pending[[ First[Ordering[angles, 1]] ]],

        rule === "loxo" && Length[pending] > 0,
          (* R2 フォールバック: max |Det| *)
          absDets = Abs[Det[{cc[[top]], cc[[last]], cc[[#]]}]] & /@ pending;
          pending[[ First[Ordering[-absDets, 1]] ]],

        rule === "maxz" && Length[pending] > 0,
          (* RZ フォールバック: min z → min Det タイブレーク *)
          zvals  = cc[[pending, 3]];
          minz   = Min[zvals];
          tied   = Select[pending, cc[[#, 3]] <= minz + 10^-10 &];
          If[Length[tied] == 1, tied[[1]],
            tDets = Det[{cc[[top]], cc[[last]], cc[[#]]}] & /@ tied;
            tied[[ First[Ordering[tDets, 1]] ]]
          ],

        True, Break[]; Null
      ];

      If[next === Null, Break[]];
      AppendTo[order, next];
      inOrder[next] = True;
      nbrs = Table[DeleteCases[nbrs[[i]], next], {i, nf}];
      prev = last;
      last = next;
    ];

    {order, Length[order] == nf}
  ];


(* ================================================================
   Section 3: 全 (top, f2) ペアの実行
   ================================================================ *)

p3lRunAll[vers_, faces_, rule_] :=
  Module[{nf, adj, lv, cc, results, f2Cands, order, complete},
    nf      = Length[faces];
    adj     = p3lBuildAdj[faces];
    results = {};
    Do[
      lv      = p3lAlignTopToZ[vers, faces, top];
      cc      = p3lFaceCentroid[#, lv, faces] & /@ Range[nf];
      f2Cands = adj[[top]];
      Do[
        {order, complete} =
          p3lPeelPair[cc, adj, nf, top, f2Cands[[f2i]], rule];
        AppendTo[results, {top, f2i, order, complete}],
        {f2i, Length[f2Cands]}
      ],
      {top, nf}
    ];
    results
  ];


(* ================================================================
   Section 4: PolyhedronData からの読み込みと処理
   ================================================================ *)

(* PolyhedronData の Graphics3D から GraphicsComplex を取り出し，
   vers (頂点座標リスト) と faces (頂点インデックスリスト) を得る．
   PolyhedronData[name,"Faces"] は Polygon オブジェクトを返すため
   直接使えない点に注意． *)
p3lExtractPolyhedron[name_] :=
  Module[{poly, gc, prim, faces},
    poly = PolyhedronData[name];
    gc   = First[Cases[poly, _GraphicsComplex, Infinity]];
    prim = gc[[2]];
    faces = Which[
      (* 全面が1つの Polygon にまとめられている場合: Polygon[{f1, f2, ...}] *)
      Head[prim] === Polygon && VectorQ[prim[[1, 1]]],
        prim[[1]],
      (* 全面が1つの Polygon + 面ごとに整数リスト *)
      Head[prim] === Polygon,
        prim[[1]],
      (* List of Polygon objects *)
      Head[prim] === List,
        Cases[prim, _Polygon, {1}][[All, 1]],
      True,
        Cases[{prim}, _Polygon, Infinity][[All, 1]]
    ];
    {N[gc[[1]]], faces}
  ];

p3lFromBuiltin[name_, rule_] :=
  Module[{vers, faces, nf, results, complete, seen, uniq},
    {vers, faces} = p3lExtractPolyhedron[name];
    nf    = Length[faces];
    Print["\n== ", name, "  (faces: ", nf, ", rule: ", rule, ") =="];
    results  = p3lRunAll[vers, faces, rule];
    complete = Select[results, #[[4]]&];
    Print["  Total (C1,C2) pairs : ", Length[results]];
    Print["  Complete             : ", Length[complete]];
    seen = <||>;
    uniq = {};
    Do[
      If[!KeyExistsQ[seen, r[[3]]],
        seen[r[[3]]] = True; AppendTo[uniq, r]],
      {r, complete}
    ];
    Print["  Unique complete      : ", Length[uniq]];
    <|"name"     -> name,
      "rule"     -> rule,
      "nFaces"   -> nf,
      "total"    -> Length[results],
      "complete" -> Length[complete],
      "unique"   -> Length[uniq],
      "details"  -> results|>
  ];


(* ================================================================
   Section 5: 全 Platonic 立体で両規則を比較
   ================================================================ *)

p3lCompareRules[] :=
  Module[{names, r2, r3, rows},
    names = {"Tetrahedron", "Cube", "Octahedron",
             "Dodecahedron", "Icosahedron"};

    Print[""];
    Print[StringRepeat["=", 55]];
    Print["  Rule 2 (loxodrome, min φ)"];
    Print[StringRepeat["=", 55]];
    r2 = Association @@ (# -> p3lFromBuiltin[#, "loxo"] & /@ names);

    Print[""];
    Print[StringRepeat["=", 55]];
    Print["  Rule 3 (max z, latitude-first)"];
    Print[StringRepeat["=", 55]];
    r3 = Association @@ (# -> p3lFromBuiltin[#, "maxz"] & /@ names);

    Print[""];
    Print[StringRepeat["=", 55]];
    Print["  比較表"];
    Print[StringRepeat["=", 55]];
    rows = Table[
      {n,
       r2[n]["nFaces"],
       r2[n]["total"],
       r2[n]["unique"], r3[n]["unique"]},
      {n, names}
    ];
    Print[TableForm[
      Prepend[rows,
        {"多面体", "面数", "総(C1,C2)",
         "ユニーク成功(規則2)", "ユニーク成功(規則3)"}]
    ]];

    <|"rule2" -> r2, "rule3" -> r3|>
  ];


Print["peeling3DLoxo.m loaded."];
Print["Usage:"];
Print["  cmp = p3lCompareRules[]             -- 全 Platonic 立体で規則2・3を比較"];
Print["  res = p3lFromBuiltin[name, rule]    -- 単独実行 (rule: \"loxo\" or \"maxz\")"];
Print["  res = p3lRunAll[vers, faces, rule]  -- 任意の多面体データで実行"];
