(* ================================================================
   gen_4d_net_figures.m
   ================================================================
   5・8・16・24胞体の代表的 3D Apple-Peel net を 2×2 グリッドで出力する．
   結果は ans4DGlobal_v4.mx から読み込む（peeling 再計算なし）．
   出力: 4d_nets.pdf
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];
Get[scriptDir <> "unfold3DExport.m"];
dataDir = scriptDir <> "_4DData/";

Print["Loading ans4DGlobal_v4.mx..."];
Get[scriptDir <> "ans4DGlobal_v4.mx"];

(* --- 色: 青 → 緑 → 赤 --- *)
apColor4D[k_, nTotal_] :=
  Blend[{RGBColor[0.20, 0.45, 0.80],
         RGBColor[0.15, 0.70, 0.40],
         RGBColor[0.90, 0.25, 0.15]},
        (k - 1) / Max[nTotal - 1, 1]];

(* --- 多胞体データ読み込み --- *)
loadPoly[n_] := Module[{raw},
  raw = Get[dataDir <> "f" <> ToString[n] <> ".m"];
  <|"vers" -> N[raw[[1]]], "faces" -> raw[[4]], "cells" -> raw[[6]]|>
];

(* 胞体ごとの視点設定 *)
viewPtFor[5]  = {1.5, -3.0,  2.5};
viewPtFor[8]  = {2.5, -2.5,  1.5};
viewPtFor[16] = {2.0, -2.0,  2.2};
viewPtFor[24] = {2.2, -1.8,  2.0};

(* --- 1多胞体の代表 net を Graphics3D として返す --- *)
makeNetGraphics[n_, labelStr_] :=
  Module[{poly, vers, faces, cells, nc, allResults, successPairs,
          localVers, order, cellLocal3Ds, cellEmbs, tris, nTotal, prims},

    poly  = loadPoly[n];
    vers  = poly["vers"];
    faces = poly["faces"];
    cells = poly["cells"];
    nc    = Length[cells];

    (* C1=1 の成功ペアから最初の1件を選ぶ *)
    allResults   = results4DGlobal[n]["results"][[1]];
    successPairs = Select[allResults, Last[#] === True &];
    If[Length[successPairs] == 0,
      Print["  WARNING: no success for C1=1, n=", n];
      Return[Graphics3D[{Text[Style["(no net)", 14], {0,0,0}]},
                        Boxed -> False, ImageSize -> {280, 260}]]
    ];

    localVers = successPairs[[1, 1]];
    order     = successPairs[[1, 2]];
    Print["  ", labelStr, ": ", Length[order], "/", nc, " cells"];

    (* 3D 展開 *)
    cellLocal3Ds = Table[
      unfoldCellLocal3D[ci, localVers, cells, faces],
      {ci, nc}
    ];
    cellEmbs = unfoldTo3DFast[cells, faces, order, cellLocal3Ds];
    tris     = unfoldTriangles[order, cells, faces, cellEmbs];
    nTotal   = Length[order];

    prims = Table[
      {FaceForm[Directive[apColor4D[k, nTotal], Opacity[0.82]]],
       EdgeForm[{Thin, GrayLevel[0.35]}],
       Polygon[tris[[k]]]},
      {k, 1, nTotal}
    ];

    Graphics3D[
      prims,
      Boxed           -> False,
      Lighting        -> "Neutral",
      SphericalRegion -> True,
      ViewPoint       -> viewPtFor[n],
      ImageSize       -> {280, 260},
      PlotLabel       -> Style[labelStr, 11, Bold]
    ]
  ];

(* ================================================================
   各胞体の図を生成
   ================================================================ *)
Print["\nGenerating nets..."];

gr5  = makeNetGraphics[5,  "5-cell"];
gr8  = makeNetGraphics[8,  "8-cell"];
gr16 = makeNetGraphics[16, "16-cell"];
gr24 = makeNetGraphics[24, "24-cell"];

(* 2×2 グリッドにまとめる *)
grAll = GraphicsGrid[
  {{gr5, gr8}, {gr16, gr24}},
  Spacings  -> {0.02, 0.02},
  ImageSize -> 600
];

outFile = scriptDir <> "4d_nets.pdf";
Export[outFile, grAll];
Print["Saved: ", outFile];
Print["Done."];
