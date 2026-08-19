(* ================================================================
   run4DExportSTLv3.m
   ================================================================
   peeling4Df4.m（v3, RZ 規則）で全5正多胞体を再計算し，
   重なり判定（checkOverlaps）を行って：
     - valid net（重なりなし）のみ STLs_v3/ に書き出す
     - valid net 比率を集計・報告する

   注：ans4DGlobal_v3.mx は 0 バイトのため，ここで直接再計算する．
   出力ディレクトリ: STLs_v3/（なければ自動作成）
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];

Get[scriptDir <> "peeling4Df4.m"];
Get[scriptDir <> "unfold3DExport.m"];

dataDir = scriptDir <> "_4DData/";
outDir  = scriptDir <> "STLs_v3/";
If[!DirectoryQ[outDir], CreateDirectory[outDir]];

polytopeNames = <|5->"5-cell", 8->"8-cell", 16->"16-cell",
                  24->"24-cell", 120->"120-cell"|>;

loadPoly[n_] := Module[{raw},
  raw = Get[dataDir <> "f" <> ToString[n] <> ".m"];
  <|"vers" -> raw[[1]], "faces" -> raw[[4]], "cells" -> raw[[6]]|>
];

(* unfoldTo3D 済みの cellEmbs を再利用して STL 出力 *)
exportSTLFromEmbs[faces_, cells_, order_, cellEmbs_, fname_] :=
  Module[{tris, gr},
    tris = unfoldTriangles[order, cells, faces, cellEmbs];
    gr   = Graphics3D[{EdgeForm[], Polygon /@ tris}];
    Export[fname, gr, "STL"];
  ];

Print[""];
Print[StringRepeat["=", 65]];
Print["  4D valid-net export  (peeling4Df4 v3 + RZ, STLs_v3/)"];
Print[StringRepeat["=", 65]];

allStats = {};

Do[
  poly  = loadPoly[n];
  vers  = N[poly["vers"]];
  faces = poly["faces"];
  cells = poly["cells"];
  nc    = Length[cells];

  Print["\n--- ", polytopeNames[n], " (", nc, " cells) ---"];

  (* ── ペーリング計算 ── *)
  t0 = AbsoluteTime[];
  ans = Table[
    If[Mod[top, 100] == 0 && nc >= 100,
      Print["  peel C1=", top, "/", nc,
            " (", Round[AbsoluteTime[] - t0, 1], " s)"]
    ];
    peeling4Df4[vers, faces, cells, top, "RZ"],
    {top, nc}
  ];
  tPeel = Round[AbsoluteTime[] - t0, 0.1];

  (* 成功 order を重複除去して収集 *)
  allOrders = DeleteDuplicates[
    Select[Flatten[ans, 1], Last[#] &][[All, 2]]
  ];
  nUniq = Length[allOrders];
  Print["  Peeling: ", tPeel, " s  |  Unique successful orders: ", nUniq];

  If[nUniq == 0,
    Print["  -> Skipping (no successful orders)."];
    AppendTo[allStats, <|"n"->n, "name"->polytopeNames[n],
                          "nUniq"->0, "nOK"->0, "nOV"->0|>];
    Continue[]
  ];

  (* ── 重なり判定 + STL 出力 ── *)
  (* 各胞の局所3D座標を事前計算（order に依存しないため 1 回のみ） *)
  cellLocal3Ds = Table[
    unfoldCellLocal3D[ci, vers, cells, faces],
    {ci, Length[cells]}
  ];
  nOK = 0; nOV = 0;
  t1  = AbsoluteTime[];

  Do[
    order    = allOrders[[i]];
    cellEmbs = unfoldTo3DFast[cells, faces, order, cellLocal3Ds];
    ovs      = checkOverlaps[order, cells, faces, cellEmbs];

    If[Length[ovs] == 0,
      (* valid net → STL 書き出し *)
      nOK++;
      fname = outDir <> "unfold_" <> ToString[n] <> "cell_" <> ToString[nOK] <> ".stl";
      exportSTLFromEmbs[faces, cells, order, cellEmbs, fname],
      nOV++
    ];

    If[Mod[i, 50] == 0 || i == nUniq,
      Print["  overlap [", i, "/", nUniq, "]  OK=", nOK, "  OV=", nOV,
            "  (", Round[AbsoluteTime[] - t1, 1], " s)"]
    ],
    {i, nUniq}
  ];

  tOver = Round[AbsoluteTime[] - t1, 0.1];
  pct   = N[nOK / nUniq * 100, 3];
  Print["  Overlap check: ", tOver, " s  |  Valid: ", nOK, "/", nUniq,
        "  (", pct, "%)  |  Overlap: ", nOV];

  AppendTo[allStats,
    <|"n"->n, "name"->polytopeNames[n],
      "nUniq"->nUniq, "nOK"->nOK, "nOV"->nOV|>],

  {n, {5, 8, 16, 24, 120}}
];

(* ── サマリ表 ─────────────────────────────────────────────── *)
Print[""];
Print[StringRepeat["=", 65]];
Print["  Summary: valid nets (no overlap) — peeling4Df4 v3 RZ"];
Print[StringRepeat["=", 65]];
Print[StringPadLeft["Polytope", 10], "  ",
      StringPadLeft["Unique",   8],  "  ",
      StringPadLeft["Valid",    7],  "  ",
      StringPadLeft["Overlap",  8],  "  ",
      StringPadLeft["Valid%",   7]];
Print[StringRepeat["-", 50]];
Do[
  s   = allStats[[i]];
  pct = If[s["nUniq"] > 0, N[s["nOK"] / s["nUniq"] * 100, 3], 0];
  Print[StringPadLeft[s["name"],          10], "  ",
        StringPadLeft[ToString[s["nUniq"]], 8], "  ",
        StringPadLeft[ToString[s["nOK"]],   7], "  ",
        StringPadLeft[ToString[s["nOV"]],   8], "  ",
        StringPadLeft[ToString[PaddedForm[pct, {5, 1}]], 7]],
  {i, Length[allStats]}
];
totalUniq  = Total[allStats[[All, "nUniq"]]];
totalValid = Total[allStats[[All, "nOK"]]];
totalPct   = If[totalUniq > 0, N[totalValid / totalUniq * 100, 3], 0];
Print[StringRepeat["-", 50]];
Print[StringPadLeft["Total", 10], "  ",
      StringPadLeft[ToString[totalUniq],  8], "  ",
      StringPadLeft[ToString[totalValid], 7], "  ",
      "  ",
      StringPadLeft[ToString[PaddedForm[totalPct, {5, 1}]], 8]];

Print["\nSTLs saved to: ", outDir];
Print["Done."];
