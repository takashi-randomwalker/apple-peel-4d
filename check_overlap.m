(* check_overlap.m
   目的: 各多胞体の全成功展開図について SAT 重なり判定を行い，
         valid net (重なりなし) の数を集計する *)

SetDirectory["/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D"];

(* 必要な関数を読み込む *)
Get["peeling4Df.m"];
Get["unfold3DExport.m"];

(* 計算済み結果を読み込む *)
Print["Loading ans4D.mx ..."];
Get["ans4D.mx"];   (* results4D を定義 *)
Print["Loaded. Keys: ", Keys[results4D]];

(* ============================================================
   各多胞体の valid net カウント
   ============================================================ *)

Print[""];
Print["==== 3D overlap check (SAT) ===="];
Print[StringRepeat["-", 60]];
Print[PaddedForm["Polytope", 12],
      PaddedForm["Orders",    9],
      PaddedForm["Valid",     8],
      PaddedForm["Overlap",   9],
      PaddedForm["Valid%",    8]];
Print[StringRepeat["-", 60]];

summaryRows = {};

Do[
  name = ToString[n] <> "-cell";
  raw   = Get["_4DData/f" <> ToString[n] <> ".m"];
  vers  = N[raw[[1]]];
  faces = raw[[4]];
  cells = raw[[6]];

  (* 成功した全連結順を収集・重複除去 *)
  allOrders = {};
  Do[
    Do[
      ord     = results4D[n]["results"][[top]][[f2]][[2]];
      success = results4D[n]["results"][[top]][[f2]][[3]];
      If[success, AppendTo[allOrders, ord]],
      {f2, Length[results4D[n]["results"][[top]]]}
    ],
    {top, Length[results4D[n]["results"]]}
  ];
  allOrders = DeleteDuplicates[allOrders];

  nOK = 0; nOV = 0;
  Do[
    cellEmbs = unfoldTo3D[vers, faces, cells, allOrders[[i]]];
    overlaps = checkOverlaps[allOrders[[i]], cells, faces, cellEmbs];
    If[Length[overlaps] == 0, nOK++, nOV++],
    {i, Length[allOrders]}
  ];

  pct = If[Length[allOrders] > 0,
    N[100 * nOK / Length[allOrders], 4], 0];

  Print[PaddedForm[name,              12],
        PaddedForm[Length[allOrders],  9],
        PaddedForm[nOK,                8],
        PaddedForm[nOV,                9],
        PaddedForm[pct,                8]];

  AppendTo[summaryRows, {name, Length[allOrders], nOK, nOV, pct}];
  , {n, Keys[results4D]}
];

Print[StringRepeat["-", 60]];
Print["Valid  : 重なりなし (valid net)"];
Print["Overlap: 重なりあり"];
