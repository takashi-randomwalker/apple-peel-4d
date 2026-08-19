(* ================================================================
   run4DExport120cellC1_3.m
   120胞体の C1=3 固定で全 C2 を計算し，成功した order を
   重なり判定なしで STLs_120cell/ に書き出す．
   出力ファイル: unfold_120cell_c1_3_1.stl, ..., _N.stl
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];

Get[scriptDir <> "peeling4Df4.m"];
Get[scriptDir <> "unfold3DExport.m"];

dataDir = scriptDir <> "_4DData/";
outDir  = scriptDir <> "STLs_120cell/";
If[!DirectoryQ[outDir], CreateDirectory[outDir]];

raw   = Get[dataDir <> "f120.m"];
vers  = N[raw[[1]]];
faces = raw[[4]];
cells = raw[[6]];
nc    = Length[cells];
Print["120-cell: ", nc, " cells, ", Length[vers], " vertices, ", Length[faces], " faces"];

fixedTop = 3;

Print["Running peeling4Df4 with C1=", fixedTop, " (RZ) ..."];
t0  = AbsoluteTime[];
ans = peeling4Df4[vers, faces, cells, fixedTop, "RZ"];
Print["  Peeling: ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* 成功した order を重複除去して収集 *)
allOrders = DeleteDuplicates[
  Select[ans, Last[#] &][[All, 2]]
];
nUniq = Length[allOrders];
Print["  Unique successful orders: ", nUniq];

If[nUniq == 0,
  Print["No successful orders. Exiting."];
  Exit[]
];

(* cellLocal3Ds を事前計算 *)
Print["Pre-computing cellLocal3Ds ..."];
t0 = AbsoluteTime[];
cellLocal3Ds = Table[unfoldCellLocal3D[ci, vers, cells, faces], {ci, nc}];
Print["  Done: ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* 全 order を STL 出力 *)
Print["Exporting ", nUniq, " STL files to ", outDir, " ..."];
t1 = AbsoluteTime[];
Do[
  order    = allOrders[[i]];
  cellEmbs = unfoldTo3DFast[cells, faces, order, cellLocal3Ds];
  tris     = unfoldTriangles[order, cells, faces, cellEmbs];
  gr       = Graphics3D[{EdgeForm[], Polygon /@ tris}];
  fname    = outDir <> "unfold_120cell_c1_3_" <> ToString[i] <> ".stl";
  Export[fname, gr, "STL"];
  If[Mod[i, 3] == 0 || i == nUniq,
    Print["  [", i, "/", nUniq, "]  (", Round[AbsoluteTime[] - t1, 1], " s)"]
  ],
  {i, nUniq}
];

Print["Done. ", nUniq, " STL files saved to: ", outDir];
