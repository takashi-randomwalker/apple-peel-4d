(* ================================================================
   gen_120cell_figures.m
   ================================================================
   Fig 1: 120-cell の 3D 展開図（C1=1, C2=2, RZ, カラーグラデーション）
   Fig 2: 120-cell の cell-centroid-up w 分布（緯度帯構造）
   出力: 120cell_unfolding.pdf, 120cell_wdist.pdf
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];
Get[scriptDir <> "unfold3DExport.m"];
Get[scriptDir <> "peeling4Df4.m"];
dataDir = scriptDir <> "_4DData/";

raw   = Get[dataDir <> "f120.m"];
vers  = N[raw[[1]]]; faces = raw[[4]]; cells = raw[[6]];
nc    = Length[cells];

Print["Computing peeling order for C1=1, C2=2..."];
t0 = AbsoluteTime[];
res = peeling4Df4[vers, faces, cells, 1, "RZ", False];
(* C2=2 は nb1 の1番目 *)
{localVers, order, ok} = res[[1]];
Print["  success=", ok, "  time=", Round[AbsoluteTime[]-t0, 0.1], " s"];

(* ============================================================
   Fig 1: 3D 展開図
   ============================================================ *)
Print["\nGenerating Fig 1: 3D unfolding..."];

cellLocal3Ds = Table[
  unfoldCellLocal3D[ci, N[vers], cells, faces],
  {ci, nc}
];
cellEmbs = unfoldTo3DFast[cells, faces, order, cellLocal3Ds];

(* 青→緑→赤グラデーション *)
apColor[k_, nTotal_] := Blend[{Blue, Green, Red}, (k-1)/(nTotal-1)];

tris = unfoldTriangles[order, cells, faces, cellEmbs];
nTotal = Length[order];
coloredPolygons = Table[
  {apColor[k, nTotal], EdgeForm[{Thin, GrayLevel[0.5]}],
   Polygon[tris[[k]]]},
  {k, 1, nTotal}
];

gr1 = Graphics3D[
  coloredPolygons,
  Boxed -> False,
  Lighting -> "Neutral",
  ImageSize -> {600, 500},
  PlotLabel -> Style["120-cell: RZ unfolding (C1=1, C2=2)", 12]
];

out1 = scriptDir <> "120cell_unfolding.pdf";
Export[out1, gr1];
Print["  Saved: ", out1];

(* ============================================================
   Fig 2: w 分布（緯度帯構造）
   ============================================================ *)
Print["\nGenerating Fig 2: w-distribution..."];

cellVerts[i_] := Union @@ faces[[ cells[[i]] ]];
wHat = {0, 0, 0, 1};

lv = (# - Mean[vers]) & /@ vers;
cc1 = Mean[lv[[ cellVerts[1] ]]];
lv = Which[
  Abs[(cc1/Norm[cc1]).wHat - 1] < 0.0001, lv,
  Abs[(cc1/Norm[cc1]).wHat + 1] < 0.0001, -lv,
  True,
    a = ArcCos[Clip[cc1.wHat/Norm[cc1], {-1,1}]];
    RotationMatrix[a, {cc1/Norm[cc1], wHat}].# & /@ lv
];
cCenters120 = Mean[lv[[ cellVerts[#] ]]] & /@ Range[nc];
wCoords = Sort[cCenters120[[All, 4]]];

(* バンドに分類 *)
eps = 0.001;
groups = {};
current = {wCoords[[1]]};
Do[
  If[Abs[wCoords[[i]] - Last[current]] < eps,
    AppendTo[current, wCoords[[i]]],
    AppendTo[groups, current]; current = {wCoords[[i]]}
  ], {i, 2, Length[wCoords]}
];
AppendTo[groups, current];

wMeans = Mean /@ groups;
counts = Length /@ groups;
nb = Length[groups];
Print["  Bands: ", nb, "  (120-cell cell-centroid-up)"];

(* 棒グラフ：w 値（横軸）× セル数（縦軸） *)
gr2 = BarChart[
  counts,
  ChartLabels -> {
    Placed[
      Table[
        If[Mod[i, 4] == 1,
          Style[ToString[Round[wMeans[[i]], 0.01]], 7],
          ""],
        {i, nb}
      ],
      Below
    ]
  },
  ChartStyle -> Table[
    Blend[{Blue, Cyan, Green, Yellow, Red},
          (wMeans[[i]] - Min[wMeans]) / (Max[wMeans] - Min[wMeans])],
    {i, nb}
  ],
  BarSpacing -> 0.1,
  Frame -> {{True, False}, {True, False}},
  FrameLabel -> {
    {Style["Number of cells", 11], None},
    {Style["Cell-centroid w-coordinate", 11], None}
  },
  PlotLabel -> Style["120-cell: cell-centroid-up latitude bands (C1=1)", 12],
  ImageSize -> {600, 300},
  GridLines -> {{}, Automatic},
  GridLinesStyle -> Directive[GrayLevel[0.8], Dashed]
];

out2 = scriptDir <> "120cell_wdist.pdf";
Export[out2, gr2];
Print["  Saved: ", out2];

(* バンド情報を表示 *)
Print["\n--- Band structure (", nb, " bands) ---"];
cumul = 0;
Do[
  cumul += counts[[i]];
  Print["  band ", StringPadLeft[ToString[i],2], ": w=",
        StringPadLeft[ToString[Round[wMeans[[i]], 0.001]], 7],
        "  n=", StringPadLeft[ToString[counts[[i]]], 3],
        "  cumul=", cumul],
  {i, nb}
];

Print["\nDone."];
