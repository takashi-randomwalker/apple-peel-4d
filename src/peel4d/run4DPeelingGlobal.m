(* ================================================================
   run4DPeelingGlobal.m
   ================================================================
   peeling4Df4（規則3 + C2 グローバル参照 xyz 行列式）を
   全6種の正多胞体に適用し，結果を Perfect / Possible / Impossible に分類する．

   比較対象:
     peeling4Df.m   (xy 2成分外積, ans4D.mx に保存済み)
     peeling4Df4.m  (C2 グローバル参照 xyz 行列式, 本スクリプト)

   使い方:
     SetDirectory["/path/to/260324Peeling4D"];
     Get["peeling4Df4.m"];
     Get["run4DPeelingGlobal.m"];

   結果は results4DGlobal に格納される．
   ans4DGlobal.mx に保存するには末尾の DumpSave 行を有効化する．
   ================================================================ *)

dataDir = If[ValueQ[$InputFileName],
             DirectoryName[$InputFileName],
             NotebookDirectory[]] <> "_4DData/";

polytopeNames = Association[
  5   -> "5-cell (pentachoron)",
  8   -> "8-cell (tesseract)",
  16  -> "16-cell (hexadecachoron)",
  24  -> "24-cell (icositetrachoron)",
  120 -> "120-cell (hecatonicosachoron)",
  600 -> "600-cell (hexacosichoron)"
];

loadPolytope4G[n_] :=
  Module[{dd, file, raw, vers, faces, cells},
    dd   = If[ValueQ[dataDir], dataDir, NotebookDirectory[] <> "_4DData/"];
    file = dd <> "f" <> ToString[n] <> ".m";
    Print["Loading: ", file];
    raw  = Get[file];
    If[!ListQ[raw] || Length[raw] != 6,
      Print["ERROR: ", file]; Return[$Failed]
    ];
    vers  = raw[[1]];
    faces = raw[[4]];
    cells = raw[[6]];
    Association["vers" -> vers, "faces" -> faces,
                "cells" -> cells, "n" -> n,
                "name" -> polytopeNames[n]]
  ];

classifyPolytope4G[data_] :=
  Module[{vers, faces, cells, n, ans, allBools,
          trueCount, falseCount, uniqueOrders, classification, t0},
    vers   = data["vers"];
    faces  = data["faces"];
    cells  = data["cells"];
    n      = data["n"];
    Print["Processing ", data["name"], " (", Length[cells], " cells) ..."];
    t0 = AbsoluteTime[];
    ans = Table[peeling4Df4[N[vers], faces, cells, top],
                {top, Length[cells]}];
    allBools   = Flatten[Map[Last, ans, {2}]];
    trueCount  = Count[allBools, True];
    falseCount = Count[allBools, False];
    uniqueOrders = DeleteDuplicates[
      Select[Flatten[ans, 1], Last[#] &][[All, 2]]
    ];
    classification = Which[
      falseCount == 0, "Perfect",
      trueCount  == 0, "Impossible",
      True,            "Possible"
    ];
    Print["  Done in ", Round[AbsoluteTime[] - t0, 0.1], " s  |  ",
          classification,
          "  (True: ", trueCount, ", False: ", falseCount,
          ", Unique: ", Length[uniqueOrders], ")"];
    Association["name" -> data["name"], "n" -> n,
                "classification" -> classification,
                "trueCount"      -> trueCount,
                "falseCount"     -> falseCount,
                "uniqueCount"    -> Length[uniqueOrders],
                "results"        -> ans]
  ];

nList = {5, 8, 16, 24, 120, 600};

Print["======================================================"];
Print["  Apple-Peel (Global C2 Ref + Rule 3) -- peeling4Df4"];
Print["======================================================"];

results4DGlobal = Association @@ Table[
  Module[{data, res},
    data = loadPolytope4G[n];
    res  = classifyPolytope4G[data];
    n -> res
  ],
  {n, nList}
];

Print[""];
Print["------ Summary ------"];
Print[TableForm[
  Prepend[
    Table[With[{r = results4DGlobal[n]},
      {r["name"], r["classification"],
       r["trueCount"], r["falseCount"], r["uniqueCount"]}],
      {n, nList}],
    {"Polytope", "Classification", "True (C1,C2)", "False (C1,C2)", "Unique"}
  ]
]];

(* 結果を保存する場合はアンコメント:
DumpSave[
  If[ValueQ[$InputFileName], DirectoryName[$InputFileName], NotebookDirectory[]]
  <> "ans4DGlobal.mx",
  results4DGlobal
]; *)

results4DGlobal
