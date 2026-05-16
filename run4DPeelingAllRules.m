(* ================================================================
   run4DPeelingAllRules.m
   ================================================================
   全6種の正多胞体に対して R1/R2/R3 の三規則を一括実行し，
   成功数・ユニーク成功数を比較する．

   前提:
     peeling4Df_rules.m が同じディレクトリにあること
     _4DData/f5.m 〜 f600.m が存在すること

   使い方 (Mathematica):
     SetDirectory["/path/to/260324Peeling4D"];
     Get["peeling4Df_rules.m"];
     Get["run4DPeelingAllRules.m"];

   出力:
     rulesResults4D  -- 全結果の Association
     "rules_results_4D.mx" に保存 (最後の行をアンコメント)
   ================================================================ *)

(* scriptDir が外部から設定済みの場合は上書きしない *)
If[!ValueQ[scriptDir],
  scriptDir = If[ValueQ[$InputFileName],
    DirectoryName[$InputFileName],
    NotebookDirectory[]]
];
dataDir = scriptDir <> "_4DData/";

(* peeling4Df_rules.m を同じディレクトリから読み込む *)
Get[scriptDir <> "peeling4Df_rules.m"];

polytopeNames = Association[
   5 -> "5-cell",
   8 -> "8-cell",
  16 -> "16-cell",
  24 -> "24-cell",
 120 -> "120-cell",
 600 -> "600-cell"
];

(* データ読み込み (run4DPeeling.m と同じ構造) *)
loadPolytope[n_] :=
  Module[{file, raw, vers, edgs, neis, faces, fneis, cells},
    file = dataDir <> "f" <> ToString[n] <> ".m";
    Print["Loading: ", file];
    raw = Get[file];
    If[!ListQ[raw] || Length[raw] != 6,
      Print["ERROR: ", file]; Return[$Failed]
    ];
    {vers, edgs, neis, faces, fneis, cells} = raw;
    Association[
      "vers"  -> vers,
      "faces" -> faces,
      "cells" -> cells,
      "n"     -> n,
      "name"  -> polytopeNames[n]
    ]
  ];

(* 1規則・1多胞体の実行 *)
runOneRule[data_, rule_] :=
  Module[{vers, faces, cells, nCells, ans, allBools,
          trueCount, falseCount, uniqueOrders, t0},
    vers   = N[data["vers"]];
    faces  = data["faces"];
    cells  = data["cells"];
    nCells = Length[cells];
    Print["  rule=", rule, "  polytope=", data["name"],
          "  (", nCells, " cells) ..."];
    t0  = AbsoluteTime[];
    ans = Table[
      peeling4DfRules[vers, faces, cells, top, rule],
      {top, nCells}
    ];
    allBools   = Flatten[Map[Last, ans, {2}]];
    trueCount  = Count[allBools, True];
    falseCount = Count[allBools, False];
    (* ユニークな成功 order を数える *)
    uniqueOrders = DeleteDuplicates[
      Select[Flatten[ans, 1], Last[#] &][[All, 2]]
    ];
    Print["    Done in ", Round[AbsoluteTime[] - t0, 0.1],
          " s  |  True: ", trueCount,
          "  Unique: ", Length[uniqueOrders]];
    Association[
      "rule"         -> rule,
      "trueCount"    -> trueCount,
      "falseCount"   -> falseCount,
      "uniqueCount"  -> Length[uniqueOrders],
      "results"      -> ans
    ]
  ];

(* 全多胞体・全規則の実行 *)
rules  = {"maxphi", "loxo", "maxw"};
(* nList が外部から設定済みの場合は上書きしない *)
If[!ValueQ[nList], nList = {5, 8, 16, 24, 120, 600}];

Print["======================================================"];
Print["  4D Apple-Peel: R1 / R2 / R3 Rule Comparison"];
Print["======================================================"];

rulesResults4D = Association @@ Table[
  Module[{data, ruleRes},
    data    = loadPolytope[n];
    ruleRes = Association @@ Table[
      r -> runOneRule[data, r],
      {r, rules}
    ];
    n -> ruleRes
  ],
  {n, nList}
];

(* サマリ表 *)
Print[""];
Print["------ Summary ------"];
Print[TableForm[
  Prepend[
    Flatten[Table[
      Table[
        With[{res = rulesResults4D[n][r]},
          {polytopeNames[n], r,
           res["trueCount"], res["uniqueCount"]}
        ],
        {r, rules}
      ],
      {n, nList}
    ], 1],
    {"Polytope", "Rule", "Success (C1,C2)", "Unique Success"}
  ]
]];

(* 結果の保存 (必要に応じてアンコメント):
DumpSave[
  If[ValueQ[$InputFileName], DirectoryName[$InputFileName], NotebookDirectory[]]
  <> "rules_results_4D.mx",
  rulesResults4D
]; *)

rulesResults4D
