(* ================================================================
   check_equivariance.m
   ================================================================
   アルキメデス13種に対して「同じ face-type ペア (F1, F2) は
   同じ成功/失敗を示すか」を検証する．

   face-type = 面の辺数（3=三角形，4=四角形，6=六角形，etc.）

   出力: 各多面体 × 各規則ごとに，全 (type(F1), type(F2)) グループが
         「全成功 / 全失敗 / 混在」のいずれかを表示．
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "peeling3DLoxo.m"];
Get[scriptDir <> "archimedean_faceup_results.mx"];

solidNames = {
  "TruncatedTetrahedron", "Cuboctahedron", "TruncatedCube",
  "TruncatedOctahedron", "Rhombicuboctahedron", "TruncatedCuboctahedron",
  "SnubCube", "Icosidodecahedron", "TruncatedDodecahedron",
  "TruncatedIcosahedron", "Rhombicosidodecahedron",
  "TruncatedIcosidodecahedron", "SnubDodecahedron"
};

rules = {"R1", "R2", "R3"};

(* 全多胞体の結果 *)
allOK = True;
summaryRows = {};

Print[""];
Print[StringRepeat["=", 70]];
Print["  Equivariance check: uniformity within (type(F1), type(F2)) groups"];
Print[StringRepeat["=", 70]];

Do[
  {vers, faces} = p3lExtractPolyhedron[name];
  (* face type = 辺数 *)
  faceType = Length /@ faces;

  solidOK = True;
  rowCols = {name};

  Do[
    details = allResults[name]["details"][rule];
    (* {top, f2, order, success} *)
    (* top, f2 は 1-indexed face番号 *)
    grouped = GroupBy[details,
      {faceType[[ #[[1]] ]], faceType[[ #[[2]] ]]} &
    ];

    ruleOK = True;
    mixedGroups = {};

    KeyValueMap[Function[{key, entries},
      successes = Count[entries, _?(#[[4]] &)];
      total     = Length[entries];
      If[successes != 0 && successes != total,
        ruleOK = False;
        AppendTo[mixedGroups,
          ToString[key] <> ": " <>
          ToString[successes] <> "/" <> ToString[total]]
      ]
    ], grouped];

    status = If[ruleOK, "OK", "MIXED"];
    If[!ruleOK,
      solidOK = False;
      allOK   = False;
      Print["  ", name, "  [", rule, "]  MIXED groups: ", mixedGroups]
    ];
    AppendTo[rowCols, status]
  , {rule, rules}];

  If[solidOK,
    Print["  ", name, "  -> all rules OK"]
  ];
  AppendTo[summaryRows, rowCols]

, {name, solidNames}];

Print[""];
Print[StringRepeat["=", 70]];
If[allOK,
  Print["RESULT: All (type(F1),type(F2)) groups are uniform for all solids and rules."];
  Print["        => No numerical-error-driven inconsistency detected."],
  Print["RESULT: Some groups are MIXED (see above)."]
];
Print[StringRepeat["=", 70]];

(* ---- 詳細: 各多面体のグループ別成功数 ---- *)
Print[""];
Print["---- Detail: success counts by (type(F1), type(F2)) ----"];
Print[""];

Do[
  {vers, faces} = p3lExtractPolyhedron[name];
  faceType = Length /@ faces;

  Print[Style[name, Bold]];
  detailsR1 = allResults[name]["details"]["R1"];
  detailsR3 = allResults[name]["details"]["R3"];

  grouped1 = GroupBy[detailsR1,
    {faceType[[ #[[1]] ]], faceType[[ #[[2]] ]]} &];
  grouped3 = GroupBy[detailsR3,
    {faceType[[ #[[1]] ]], faceType[[ #[[2]] ]]} &];

  allKeys = Union[Keys[grouped1], Keys[grouped3]];
  Do[
    n1 = Length[grouped1[key]];
    s1 = Count[grouped1[key], _?(#[[4]] &)];
    n3 = If[KeyExistsQ[grouped3, key], Length[grouped3[key]], 0];
    s3 = If[KeyExistsQ[grouped3, key], Count[grouped3[key], _?(#[[4]] &)], 0];
    Print["  type(F1,F2)=", key,
          "  n=", n1,
          "  RS: ", s1, "/", n1,
          "  RZ: ", s3, "/", n3]
  , {key, allKeys}];
  Print[""]

, {name, solidNames}];
