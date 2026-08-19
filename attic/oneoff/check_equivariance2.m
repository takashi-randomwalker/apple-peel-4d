(* ================================================================
   check_equivariance2.m
   ================================================================
   allResults の {top, f2i, order, success} において
   f2i = adj[[top]] 内の位置番号であることを考慮し，
   order[[1]], order[[2]] から実際の面インデックスを取得して
   face-type グループを正しく求める．
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

allOK   = True;
detailTable = {};

Print[""];
Print[StringRepeat["=", 70]];
Print["  Equivariance check (corrected): face-type group uniformity"];
Print[StringRepeat["=", 70]];

Do[
  {vers, faces} = p3lExtractPolyhedron[name];
  faceType = Length /@ faces;
  solidOK  = True;

  Do[
    details = allResults[name]["details"][rule];
    (* order[[1]] = F1 面インデックス（= top と一致するはず）       *)
    (* order[[2]] = F2 実面インデックス                             *)
    (* order が長さ1以下のペアは除外（理論上ありえないが安全のため） *)
    validDetails = Select[details, Length[#[[3]]] >= 2 &];

    grouped = GroupBy[validDetails,
      {faceType[[ #[[3, 1]] ]], faceType[[ #[[3, 2]] ]]} &
    ];

    ruleOK     = True;
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

    If[!ruleOK,
      solidOK = False; allOK = False;
      Print["  [MIXED] ", name, "  [", rule, "] ", mixedGroups]
    ]
  , {rule, rules}];

  If[solidOK, Print["  [OK]   ", name]]

, {name, solidNames}];

Print[""];
Print[StringRepeat["=", 70]];
If[allOK,
  Print["RESULT: All face-type groups are uniform. No numerical inconsistency."],
  Print["RESULT: Some MIXED groups remain (see above)."]
];
Print[StringRepeat["=", 70]];

(* ---- 詳細表: 各多面体 × (type(F1), type(F2)) グループの成功数 ---- *)
Print[""];
Print["---- Detail: RS and RZ success by (type(F1), type(F2)) ----"];
Print[""];

Do[
  {vers, faces} = p3lExtractPolyhedron[name];
  faceType = Length /@ faces;

  detR1 = allResults[name]["details"]["R1"];
  detR3 = allResults[name]["details"]["R3"];

  validR1 = Select[detR1, Length[#[[3]]] >= 2 &];
  validR3 = Select[detR3, Length[#[[3]]] >= 2 &];

  gr1 = GroupBy[validR1, {faceType[[#[[3,1]]]], faceType[[#[[3,2]]]]} &];
  gr3 = GroupBy[validR3, {faceType[[#[[3,1]]]], faceType[[#[[3,2]]]]} &];

  allKeys = Union[Keys[gr1], Keys[gr3]];

  Print[Style[name, Bold], "  (face types: ", Counts[faceType], ")"];
  Do[
    n1 = If[KeyExistsQ[gr1, k], Length[gr1[k]], 0];
    s1 = If[KeyExistsQ[gr1, k], Count[gr1[k], _?(#[[4]]&)], 0];
    n3 = If[KeyExistsQ[gr3, k], Length[gr3[k]], 0];
    s3 = If[KeyExistsQ[gr3, k], Count[gr3[k], _?(#[[4]]&)], 0];
    Print["  (F1-type, F2-type)=", k,
          "  n=", n1,
          "  RS: ", s1, "/", n1,
          "  RZ: ", s3, "/", n3]
  , {k, allKeys}];
  Print[]
, {name, solidNames}];
