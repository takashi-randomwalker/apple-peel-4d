(* analyze_120cell_c2types.m
   C1=3 固定で 120胞体の C2 候補を face-up 後 w 座標で分類 *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];

Get[scriptDir <> "peeling4Df4.m"];

dataDir = scriptDir <> "_4DData/";
raw   = Get[dataDir <> "f120.m"];
vers  = N[raw[[1]]];
faces = raw[[4]];
cells = raw[[6]];
nc    = Length[cells];

fixedTop = 3;
wHat     = {0., 0., 0., 1.};

(* 全セルの頂点インデックス *)
cellVerts = Function[ci,
  Union[Flatten[Map[faces[[#]] &, cells[[ci]]]]]
];

(* 元座標での全セル重心 *)
cc0 = Table[Mean[vers[[cellVerts[ci]]]], {ci, nc}];

(* face-up 回転: C1 重心方向 → +w 軸 *)
cc1 = cc0[[fixedTop]];
rotMat =
  Which[
    Abs[(cc1 / Norm[cc1]) . wHat - 1] < 0.0001, IdentityMatrix[4],
    Abs[(cc1 / Norm[cc1]) . wHat + 1] < 0.0001, -IdentityMatrix[4],
    True,
      RotationMatrix[ArcCos[Clip[cc1 . wHat / Norm[cc1], {-1,1}]],
                     {cc1 / Norm[cc1], wHat}]
  ];

(* face-up 後の全セル重心 *)
cc = Map[rotMat . # &, cc0];

(* C2 候補: C1 と面を共有するセル *)
c2cands = Select[Range[nc],
  # != fixedTop &&
  Length[Intersection[cells[[fixedTop]], cells[[#]]]] >= 1 &
];
Print["C1=", fixedTop, "  C2 candidates: ", Length[c2cands]];

(* w 座標でソート *)
c2w = SortBy[Map[{#, cc[[#, 4]]} &, c2cands], -Last[#] &];
Print["\nC2 candidates by face-up w (high → low):"];
Do[
  Print["  C2=", c2w[[i,1]], "  w=", NumberForm[c2w[[i,2]], {8,5}]],
  {i, Length[c2w]}
];

(* w 値でグループ化（差 < 10^-4 を同グループ） *)
wVals  = c2w[[All, 2]];
groups = {};
used   = ConstantArray[False, Length[wVals]];
Do[
  If[!used[[i]],
    grp = {c2w[[i, 1]]};
    used[[i]] = True;
    Do[
      If[!used[[j]] && Abs[wVals[[j]] - wVals[[i]]] < 10^-4,
        AppendTo[grp, c2w[[j, 1]]];
        used[[j]] = True
      ],
      {j, i+1, Length[wVals]}
    ];
    AppendTo[groups, {wVals[[i]], grp}]
  ],
  {i, Length[wVals]}
];

Print["\nLatitude bands:"];
Do[
  g = groups[[k]];
  Print["  w = ", NumberForm[g[[1]], {8,5}],
        "  (", Length[g[[2]]], " cells): C2 = ", g[[2]]],
  {k, Length[groups]}
];

(* peeling 実行 *)
Print["\nRunning peeling for C1=", fixedTop, " (RZ) ..."];
t0  = AbsoluteTime[];
ans = peeling4Df4[vers, faces, cells, fixedTop, "RZ"];
Print["Peeling: ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* 成功/失敗ごとに C2 = order[[2]] を取得 *)
successC2 = Select[ans, Last[#] &][[All, 2, 2]];  (* order[[2]] = C2 index *)
failC2    = Select[ans, !Last[#] &][[All, 2, 2]];
Print["Success C2 (", Length[successC2], "): ", Sort[successC2]];
Print["Fail    C2 (", Length[failC2],    "): ", Sort[failC2]];

(* グループ × 成功/失敗 *)
Print["\nBand × success:"];
Do[
  g = groups[[k]];
  ok = Select[g[[2]], MemberQ[successC2, #] &];
  ng = Select[g[[2]], MemberQ[failC2,    #] &];
  Print["  w = ", NumberForm[g[[1]], {8,5}],
        "  size=", Length[g[[2]]],
        "  OK=", Sort[ok],
        "  FAIL=", Sort[ng]],
  {k, Length[groups]}
];

(* STL ファイルとの対応（i番目の成功order → unfold_120cell_c1_3_i.stl） *)
orderedSucc = Select[ans, Last[#] &];
orderedSucc = DeleteDuplicates[orderedSucc, #1[[2]] === #2[[2]] &];
Print["\nSTL file → C2 mapping:"];
Do[
  c2 = orderedSucc[[i, 2, 2]];
  wc2 = cc[[c2, 4]];
  Print["  unfold_120cell_c1_3_", i, ".stl  C2=", c2,
        "  w=", NumberForm[wc2, {8,5}]],
  {i, Length[orderedSucc]}
];

Print["Done."];
