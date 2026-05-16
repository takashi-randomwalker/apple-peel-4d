(* ================================================================
   run4DExportSTLv4.m
   ================================================================
   ans4DGlobal_v4.mx（k>=3 Det=0 フォールバック修正版）から
   成功 order を読み込み，重なり判定後 valid net を STLs_v4/ に出力．

   v3 との違い:
     - ペーリング再計算なし（ans4DGlobal_v4.mx から order を取得）
     - 120胞体 Perfect 1440 orders に対応
   出力ディレクトリ: STLs_v4/
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];

Get[scriptDir <> "unfold3DExport.m"];

dataDir = scriptDir <> "_4DData/";
outDir  = scriptDir <> "STLs_v4/";
If[!DirectoryQ[outDir], CreateDirectory[outDir]];

(* 計算結果を読み込む *)
v4File = scriptDir <> "ans4DGlobal_v4.mx";
If[FileExistsQ[v4File],
  Get[v4File]; Print["Loaded ans4DGlobal_v4.mx"],
  Print["ERROR: ans4DGlobal_v4.mx not found."]; Abort[]
];

polytopeNames = <|5->"5-cell", 8->"8-cell", 16->"16-cell",
                  24->"24-cell", 120->"120-cell"|>;

loadPoly[n_] := Module[{raw},
  raw = Get[dataDir <> "f" <> ToString[n] <> ".m"];
  <|"vers" -> raw[[1]], "faces" -> raw[[4]], "cells" -> raw[[6]]|>
];

exportSTLFromEmbs[faces_, cells_, order_, cellEmbs_, fname_] :=
  Module[{tris, gr},
    tris = unfoldTriangles[order, cells, faces, cellEmbs];
    gr   = Graphics3D[{EdgeForm[], Polygon /@ tris}];
    Export[fname, gr, "STL"];
  ];

Print[""];
Print[StringRepeat["=", 65]];
Print["  4D valid-net export  (peeling4Df4 v4, STLs_v4/)"];
Print[StringRepeat["=", 65]];

allStats = {};

Do[
  If[!KeyExistsQ[results4DGlobal, n],
    Print["\n--- ", polytopeNames[n], ": no data, skipping ---"];
    Continue[]
  ];

  poly  = loadPoly[n];
  faces = poly["faces"];
  cells = poly["cells"];
  nc    = Length[cells];

  Print["\n--- ", polytopeNames[n], " (", nc, " cells) ---"];

  (* 成功 order を v4 結果から取得（再計算不要） *)
  allOrders = DeleteDuplicates[
    Select[Flatten[results4DGlobal[n]["results"], 1], Last[#] &][[All, 2]]
  ];
  nUniq = Length[allOrders];
  Print["  Unique successful orders: ", nUniq,
        "  (classification: ", results4DGlobal[n]["classification"], ")"];

  If[nUniq == 0,
    Print["  -> Skipping (no successful orders)."];
    AppendTo[allStats, <|"n"->n, "name"->polytopeNames[n],
                          "nUniq"->0, "nOK"->0, "nOV"->0|>];
    Continue[]
  ];

  (* 各セルの局所3D座標を事前計算 *)
  cellLocal3Ds = Table[
    unfoldCellLocal3D[ci, N[poly["vers"]], cells, faces],
    {ci, Length[cells]}
  ];

  nOK = 0; nOV = 0;
  t1  = AbsoluteTime[];

  Do[
    order    = allOrders[[i]];
    cellEmbs = unfoldTo3DFast[cells, faces, order, cellLocal3Ds];
    ovs      = checkOverlaps[order, cells, faces, cellEmbs];

    If[Length[ovs] == 0,
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
  pct   = If[nUniq > 0, N[nOK / nUniq * 100, 3], 0];
  Print["  Overlap check: ", tOver, " s  |  Valid: ", nOK, "/", nUniq,
        "  (", pct, "%)"];

  AppendTo[allStats,
    <|"n"->n, "name"->polytopeNames[n],
      "nUniq"->nUniq, "nOK"->nOK, "nOV"->nOV|>],

  {n, {5, 8, 16, 24, 120}}
];

(* サマリ表 *)
Print[""];
Print[StringRepeat["=", 65]];
Print["  Summary: valid nets (no overlap) — peeling4Df4 v4 RZ"];
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
