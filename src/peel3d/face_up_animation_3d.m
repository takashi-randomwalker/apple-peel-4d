(* ================================================================
   face_up_animation_3d.m
   ----------------------------------------------------------------
   3次元多面体の「face-up 回転 → Apple-Peel 展開」アニメーション．

   Phase 1: 開始面 F1 の重心を +z 軸に揃える回転（= face-up）
   Phase 2: faceOrder に従い，各ヒンジ辺（hingeEdgeIds）周りに
            二面角ぶん順に折り倒して平面ネットへ展開

   依存: face_up_sequence_3d.m -> peeling3DLoxo.m（同ディレクトリ）

   使い方（末尾の driver がそのまま実行）:
     Get["face_up_animation_3d.m"];
   個別に呼ぶ場合:
     frames = p3lBuildFrames["Dodecahedron", 1, "maxz"];
     Export["out.gif", frames, "DisplayDurations"->0.05];
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
If[scriptDir === "", scriptDir = NotebookDirectory[]];
Get[scriptDir <> "face_up_sequence_3d.m"];


(* ----------------------------------------------------------------
   設定パラメータ（必要に応じて変更）
   ---------------------------------------------------------------- *)
p3lAnimName     = "Dodecahedron";  (* 対象多面体 *)
p3lAnimTop      = 1;               (* 開始面 F1 *)
p3lAnimRule     = "maxz";          (* "maxz"(RZ) / "maxphi"(RS) / "loxo"(R2) *)
p3lAnimNPhase1  = 24;              (* face-up 回転のフレーム数 *)
p3lAnimNPerFold = 12;              (* 1回の折りのフレーム数 *)
p3lAnimPause    = 8;               (* 各フェーズ端での静止フレーム数 *)
p3lAnimSize     = 500;             (* 画像サイズ (px) *)
p3lAnimView     = {1.6, -2.6, 1.9};(* 視点 *)
p3lAnimFPS      = 20;              (* 出力フレームレート *)


(* 色設定 *)
p3lColBase    = RGBColor[0.20, 0.45, 0.85];  (* F1（開始面）: 青 *)
p3lColDone    = RGBColor[0.90, 0.78, 0.15];  (* 展開済み: 黄 *)
p3lColActive  = RGBColor[0.95, 0.45, 0.10];  (* 折り曲げ中: 橙 *)
p3lColPending = RGBColor[0.25, 0.70, 0.30];  (* 未展開: 緑 *)


(* ----------------------------------------------------------------
   面の外向き法線（符号は固定 signFix で保持）
   ---------------------------------------------------------------- *)
p3lRawNormal[poly_] :=
  Normalize[Cross[poly[[2]] - poly[[1]], poly[[3]] - poly[[1]]]];


(* ----------------------------------------------------------------
   展開状態の計算
     verFU  : 各面の頂点座標（face-up 後，nf 個のポリゴン）
     order  : 面の展開順
     faces  : 面→頂点番号
     signFix: 各面の法線符号（外向き正）
     s      : 展開パラメータ 0（完全に立体）〜 (n-1)（完全に平面）
   戻り値: 各面の現在座標（nf 個のポリゴン）
   ---------------------------------------------------------------- *)
p3lUnfoldState[verFU_, order_, faces_, signFix_, s_] :=
  Module[{pos, n, K, frac, coordOf, doFold},
    pos = verFU;
    n   = Length[order];
    K   = Floor[s];
    frac = s - K;

    (* 面 f 内の大域頂点 g の現在座標 *)
    coordOf[f_, g_] := pos[[f]][[ First@First@Position[faces[[f]], g] ]];

    (* fold m: order[[m+1..n]] を hinge(order[[m]],order[[m+1]]) 周りに frac*角度 回す *)
    doFold[m_, f_] := Module[{a, b, sh, p1, p2, axis, pivot, nb, nc,
                              cosA, sinA, ang, rot, t},
      a = order[[m]]; b = order[[m + 1]];
      sh = Intersection[faces[[a]], faces[[b]]];
      If[Length[sh] < 2, Return[]];
      sh = sh[[1 ;; 2]];
      p1 = coordOf[a, sh[[1]]];   (* ヒンジ端点（base 側 = 不動側から取得） *)
      p2 = coordOf[a, sh[[2]]];
      axis  = Normalize[p2 - p1];
      pivot = p1;
      nb = signFix[[a]] p3lRawNormal[pos[[a]]];   (* base 面の外向き法線 *)
      nc = signFix[[b]] p3lRawNormal[pos[[b]]];   (* child 面の外向き法線 *)
      cosA = Clip[nc . nb, {-1., 1.}];
      sinA = axis . Cross[nc, nb];
      ang  = ArcTan[cosA, sinA];   (* nc を nb に重ねる符号付き角 *)
      rot  = RotationTransform[f ang, axis, pivot];
      Do[pos[[order[[t]]]] = rot /@ pos[[order[[t]]]], {t, m + 1, n}];
    ];

    (* 完了した折り *)
    Do[doFold[m, 1.], {m, 1, K}];
    (* 進行中の折り *)
    If[K + 1 <= n - 1 && frac > 0, doFold[K + 1, frac]];
    pos
  ];


(* ----------------------------------------------------------------
   共通描画枠（視点・範囲固定）
   ---------------------------------------------------------------- *)
p3lFrame[prims_, plotRange_] :=
  Graphics3D[prims,
    ViewPoint    -> p3lAnimView,
    ViewVertical -> {0, 0, 1},
    PlotRange    -> plotRange,
    Boxed        -> False,
    Lighting     -> "Neutral",
    SphericalRegion -> True,
    ImageSize    -> p3lAnimSize,
    Background   -> White];


(* ----------------------------------------------------------------
   全フレーム生成
   ---------------------------------------------------------------- *)
p3lBuildFrames[name_, top_, rule_, f2_: Automatic] :=
  Module[{res, faces, order, nf, versC, R, verFUsolid, signFix,
          logR, RtOf, edgeSty, drawFaces, drawSolid, frames,
          finalNet, allPts, pr, coordOf, hingeLine, sMax, mIdx,
          zAxisPrim, topDir},

    (* --- シーケンス情報 --- *)
    res    = p3lPeelSequence[name, top, rule, f2];
    faces  = res["faces"];
    order  = res["faceOrder"];
    nf     = res["nFaces"];
    versC  = res["verticesCentered"];
    topDir = res["topCentroidDir"];

    (* face-up 回転（単一の proper rotation で phase1 をスラープ） *)
    R = res["faceUpRotation"];
    logR = Re[MatrixLog[R]];
    RtOf[t_] := Re[MatrixExp[t logR]];

    (* face-up 後の各面ポリゴン（phase2 の基準座標） *)
    verFUsolid = Map[(versC . Transpose[R])[[#]] &, faces];

    (* 法線符号（原点=立体中心，外向き=法線・重心>0） *)
    signFix = Table[
      With[{poly = verFUsolid[[i]]},
        If[p3lRawNormal[poly] . Mean[poly] < 0, -1, 1]],
      {i, nf}];

    edgeSty = EdgeForm[{GrayLevel[0.25], AbsoluteThickness[1.1]}];

    (* --- 立体（phase1）を任意回転 t で描く --- *)
    drawSolid[t_] := Module[{coords, vf, arrowTo},
      coords = versC . Transpose[RtOf[t]];
      vf = Map[coords[[#]] &, faces];
      arrowTo = 1.7 (topDir . Transpose[RtOf[t]]);
      {
        Table[{If[i == order[[1]], p3lColBase, p3lColPending],
               edgeSty, Polygon[vf[[i]]]}, {i, nf}],
        {Arrowheads[0.04], AbsoluteThickness[2.2], Black,
         Arrow[{{0, 0, 0}, arrowTo}]}
      }];

    (* --- 展開状態 s を描く（色分け・現在ヒンジ強調） --- *)
    sMax = nf - 1;
    drawFaces[s_] := Module[{pos, K, colFor, prims, a, b, sh, p1, p2},
      pos = p3lUnfoldState[verFUsolid, order, faces, signFix, s];
      K = Floor[s];
      colFor[j_] := Which[
        j == 1,               p3lColBase,
        j <= K + 1,           p3lColDone,
        j == K + 2 && s > K,  p3lColActive,
        True,                 p3lColPending];
      prims = Table[
        {colFor[j], edgeSty, Polygon[pos[[ order[[j]] ]]]},
        {j, nf}];
      (* 現在のヒンジ辺を赤で強調 *)
      If[K + 1 <= sMax,
        a = order[[K + 1]]; b = order[[K + 2]];
        sh = Intersection[faces[[a]], faces[[b]]][[1 ;; 2]];
        p1 = pos[[a]][[ First@First@Position[faces[[a]], sh[[1]]] ]];
        p2 = pos[[a]][[ First@First@Position[faces[[a]], sh[[2]]] ]];
        prims = Append[prims,
          {RGBColor[0.85, 0.1, 0.1], AbsoluteThickness[4], Line[{p1, p2}]}]
      ];
      prims];

    (* --- PlotRange を最終ネットと立体の和で固定 --- *)
    finalNet = p3lUnfoldState[verFUsolid, order, faces, signFix, N[sMax]];
    allPts = Join[Flatten[verFUsolid, 1], Flatten[finalNet, 1]];
    pr = Transpose[{Min /@ Transpose[allPts], Max /@ Transpose[allPts]}];
    pr = MapThread[{#1 - 0.12 (#2 - #1), #2 + 0.12 (#2 - #1)} &,
                   {pr[[All, 1]], pr[[All, 2]]}];

    (* --- フレーム列 --- *)
    frames = {};
    (* 静止（立体） *)
    Do[AppendTo[frames, p3lFrame[drawSolid[0.], pr]], {p3lAnimPause}];
    (* face-up 回転 *)
    Do[AppendTo[frames, p3lFrame[drawSolid[t], pr]],
       {t, Subdivide[0., 1., p3lAnimNPhase1]}];
    Do[AppendTo[frames, p3lFrame[drawSolid[1.], pr]], {p3lAnimPause}];
    (* 展開（face-up 後の立体 s=0 から） *)
    Do[AppendTo[frames, p3lFrame[drawFaces[0.], pr]], {p3lAnimPause}];
    Do[AppendTo[frames, p3lFrame[drawFaces[s], pr]],
       {s, Subdivide[0., N[sMax], sMax*p3lAnimNPerFold]}];
    Do[AppendTo[frames, p3lFrame[drawFaces[N[sMax]], pr]], {2 p3lAnimPause}];

    Print["  frames: ", Length[frames],
          "  (faces=", nf, ", complete=", res["complete"], ")"];
    frames
  ];


(* ----------------------------------------------------------------
   Driver: GIF（＋可能なら MOV）を書き出す
   ---------------------------------------------------------------- *)
p3lExportAnimation[name_, top_, rule_, f2_: Automatic] :=
  Module[{frames, base, gif, dur},
    Print["Building animation: ", name, " (top=", top, ", rule=", rule, ")"];
    frames = p3lBuildFrames[name, top, rule, f2];
    dur = 1./p3lAnimFPS;
    base = scriptDir <> "faceup_" <> name <> "_" <> rule;
    gif = base <> ".gif";
    Export[gif, frames, "DisplayDurations" -> dur];
    Print["Exported GIF: ", gif];
    (* MOV は環境依存。失敗しても GIF は出力済み *)
    Quiet@Check[
      Export[base <> ".mov", frames, "FrameRate" -> p3lAnimFPS];
      Print["Exported MOV: ", base <> ".mov"],
      Print["MOV export not available (GIF only). ffmpeg で変換可: ",
            "ffmpeg -i ", gif, " ", base, ".mov"]
    ];
    gif
  ];


(* ================================================================
   実行
   ================================================================ *)
p3lExportAnimation[p3lAnimName, p3lAnimTop, p3lAnimRule];
