(* ================================================================
   run4DPeelingR3.m
   ================================================================
   peeling4Df3（規則3 + 行列式左条件）を全6正多胞体に適用し，
   結果を Perfect / Possible / Impossible に分類する．

   使い方:
     SetDirectory["/path/to/260324Peeling4D"];
     Get["peeling4Df3.m"];
     Get["run4DPeelingR3.m"];

   結果は results4DR3 に格納される．
   ans4DR3.mx に保存するには最後の DumpSave 行を有効化する．
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

loadPolytope[n_] :=
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

classifyPolytopeR3[data_] :=
  Module[{vers, faces, cells, n, ans, allBools,
          trueCount, falseCount, classification, t0},
    vers   = data["vers"];
    faces  = data["faces"];
    cells  = data["cells"];
    n      = data["n"];
    Print["Processing ", data["name"], " (", Length[cells], " cells) ..."];
    t0 = AbsoluteTime[];
    ans = Table[peeling4Df3[N[vers], faces, cells, top],
                {top, Length[cells]}];
    allBools   = Flatten[Map[Last, ans, {2}]];
    trueCount  = Count[allBools, True];
    falseCount = Count[allBools, False];
    classification = Which[
      falseCount == 0, "Perfect",
      trueCount  == 0, "Impossible",
      True,            "Possible"
    ];
    Print["  Done in ", Round[AbsoluteTime[] - t0, 0.1], " s  |  ",
          classification,
          "  (True: ", trueCount, ", False: ", falseCount, ")"];
    Association["name" -> data["name"], "n" -> n,
                "classification" -> classification,
                "trueCount" -> trueCount, "falseCount" -> falseCount,
                "results" -> ans]
  ];

nList = {5, 8, 16, 24, 120, 600};

Print["======================================================"];
Print["  Apple-Peel (Rule 3 + Det left) -- peeling4Df3"];
Print["======================================================"];

results4DR3 = Association @@ Table[
  Module[{data, res},
    data = loadPolytope[n];
    res  = classifyPolytopeR3[data];
    n -> res
  ],
  {n, nList}
];

Print["------ Summary ------"];
Print[TableForm[
  Prepend[
    Table[With[{r = results4DR3[n]},
      {r["name"], r["classification"], r["trueCount"], r["falseCount"]}],
      {n, nList}],
    {"Polytope", "Classification", "True (C1,C2)", "False (C1,C2)"}
  ]
]];

(* Uncomment to save:
DumpSave[NotebookDirectory[] <> "ans4DR3.mx", results4DR3]; *)

results4DR3
