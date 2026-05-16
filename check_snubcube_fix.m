(* ================================================================
   check_snubcube_fix.m
   ================================================================
   peeling3DLoxo.m の antipodal-case 修正を SnubCube で検証する．
   archimedean_faceup_results.mx の SnubCube エントリを更新して
   equivariance check を実行する．
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "peeling3DLoxo.m"];
Get[scriptDir <> "archimedean_faceup_results.mx"];

name = "SnubCube";
{vers, faces} = p3lExtractPolyhedron[name];
nf    = Length[faces];
adj   = p3lBuildAdj[faces];
totalPairs = Total[Length /@ adj];

Print["Re-computing ", name, " with fixed p3lAlignTopToZ ..."];

r1 = p3lRunAll[vers, faces, "maxphi"];
r2 = p3lRunAll[vers, faces, "loxo"];
r3 = p3lRunAll[vers, faces, "maxz"];
n1 = Length[Select[r1, #[[4]]&]];
n2 = Length[Select[r2, #[[4]]&]];
n3 = Length[Select[r3, #[[4]]&]];

Print["  R1 (RS): ", n1, "/", totalPairs, "  (", NumberForm[N[100 n1/totalPairs], {4,1}], "%)"];
Print["  R2:      ", n2, "/", totalPairs, "  (", NumberForm[N[100 n2/totalPairs], {4,1}], "%)"];
Print["  R3 (RZ): ", n3, "/", totalPairs, "  (", NumberForm[N[100 n3/totalPairs], {4,1}], "%)"];

(* archimedean_faceup_results.mx を読み込んだ allResults を更新 *)
allResults[name] = <|
  "nFaces" -> nf, "totalPairs" -> totalPairs,
  "R1" -> n1, "R2" -> n2, "R3" -> n3,
  "details" -> <|"R1" -> r1, "R2" -> r2, "R3" -> r3|>
|>;

(* --- Triangle subtype analysis --- *)
faceType      = Length /@ faces;
squareFaces   = Select[Range[nf], faceType[[#]] == 4 &];
triFaces      = Select[Range[nf], faceType[[#]] == 3 &];
triAdjSqCount = Association @ Table[
  f -> Length[Intersection[adj[[f]], squareFaces]],
  {f, triFaces}
];
snubTri   = Select[triFaces, triAdjSqCount[#] == 1 &];
gyrateTri = Select[triFaces, triAdjSqCount[#] == 0 &];
Print["\nTriangle subtypes: snub=", Length[snubTri], "  gyrate=", Length[gyrateTri]];

getSubtype[f_] := Which[
  faceType[[f]] == 4,      "sq",
  MemberQ[snubTri, f],     "snub",
  MemberQ[gyrateTri, f],   "gyrate"
];

(* --- Sub-orbit breakdown after fix --- *)
Print["\n--- Sub-orbit breakdown RS (after fix) ---"];
pairsR1 = Select[r1, Length[#[[3]]] >= 2 &];
pairWithSubtype = {getSubtype[#[[3,1]]], getSubtype[#[[3,2]]], #[[4]]} & /@ pairsR1;
grouped = GroupBy[pairWithSubtype, {#[[1]], #[[2]]} &];
Do[
  entries = grouped[k];
  nsuc = Count[entries, _?(#[[3]] &)];
  ntot = Length[entries];
  Print["  ", k, "  n=", ntot, "  success: ", nsuc, "/", ntot]
, {k, Keys[grouped]}];

Print["\n--- Sub-orbit breakdown RZ (after fix) ---"];
pairsR3 = Select[r3, Length[#[[3]]] >= 2 &];
pairWithSubtypeR3 = {getSubtype[#[[3,1]]], getSubtype[#[[3,2]]], #[[4]]} & /@ pairsR3;
groupedR3 = GroupBy[pairWithSubtypeR3, {#[[1]], #[[2]]} &];
Do[
  entries = groupedR3[k];
  nsuc = Count[entries, _?(#[[3]] &)];
  ntot = Length[entries];
  Print["  ", k, "  n=", ntot, "  success: ", nsuc, "/", ntot]
, {k, Keys[groupedR3]}];

(* --- Face-type equivariance check for SnubCube --- *)
Print["\n--- Face-type equivariance check (SnubCube) ---"];
rules = {"R1", "R2", "R3"};
allOK = True;
Do[
  details = allResults[name]["details"][rule];
  validDetails = Select[details, Length[#[[3]]] >= 2 &];
  grouped2 = GroupBy[validDetails,
    {faceType[[ #[[3, 1]] ]], faceType[[ #[[3, 2]] ]]} &
  ];
  ruleOK = True;
  KeyValueMap[Function[{key, entries},
    successes = Count[entries, _?(#[[4]] &)];
    total     = Length[entries];
    If[successes != 0 && successes != total,
      ruleOK = False; allOK = False;
      Print["  [MIXED] ", rule, "  ", key, ": ", successes, "/", total]
    ]
  ], grouped2];
  If[ruleOK, Print["  [OK]   ", rule]]
, {rule, rules}];

If[allOK,
  Print["\nRESULT: SnubCube face-type groups are uniform. Fix validated."],
  Print["\nRESULT: Still MIXED groups after fix — further investigation needed."]
];
