(* ================================================================
   check_snubdodec.m
   ================================================================
   SnubDodecahedron のサブ軌道分析と antipodal-case の有無を調べる．
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "peeling3DLoxo.m"];
Get[scriptDir <> "archimedean_faceup_results.mx"];

name = "SnubDodecahedron";
{vers, faces} = p3lExtractPolyhedron[name];
nf    = Length[faces];
adj   = p3lBuildAdj[faces];
totalPairs = Total[Length /@ adj];
faceType = Length /@ faces;

Print["\n=== ", name, " ==="];
Print["Faces: ", nf, "  face-type counts: ", Normal[Counts[faceType]]];
Print["Total pairs: ", totalPairs];

(* --- Triangle subtypes --- *)
pentFaces = Select[Range[nf], faceType[[#]] == 5 &];
triFaces  = Select[Range[nf], faceType[[#]] == 3 &];

triAdjPentCount = Association @ Table[
  f -> Length[Intersection[adj[[f]], pentFaces]],
  {f, triFaces}
];
triAdjCounts = Values[triAdjPentCount];
Print["\nTriangle adj-to-pentagon count distribution: ", Tally[Sort[triAdjCounts]]];

snubTri   = Select[triFaces, triAdjPentCount[#] >= 1 &];
gyrateTri = Select[triFaces, triAdjPentCount[#] == 0 &];
Print["snub triangles (adj>=1 pentagon): ", Length[snubTri]];
Print["gyrate triangles (adj=0 pentagon): ", Length[gyrateTri]];

getSubtype[f_] := Which[
  faceType[[f]] == 5,      "pent",
  MemberQ[snubTri, f],     "snub",
  MemberQ[gyrateTri, f],   "gyrate"
];

(* --- Antipodal-case 検出 --- *)
(* p3lAlignTopToZ の antipodal 条件: top重心 dot (0,0,1) <= -1+1e-4 *)
Print["\n--- Antipodal-case detection ---"];
centeredVers = N[vers - ConstantArray[Mean[vers], Length[vers]]];
faceCentroids = Mean[centeredVers[[ # ]]] & /@ faces;
dots = Normalize[#] . {0., 0., 1.} & /@ faceCentroids;
antipodalFaces = Select[Range[nf], Abs[dots[[#]] + 1] <= 10^-4 &];
Print["Faces with centroid near -z (antipodal trigger): ", antipodalFaces];
If[Length[antipodalFaces] > 0,
  Do[
    Print["  face ", f, "  type=", faceType[[f]], "  dot=", dots[[f]],
          "  centroid=", N[faceCentroids[[f]], 4]]
  , {f, antipodalFaces}],
  Print["  None found — antipodal case not triggered in any top-face alignment"]
];

(* --- sub-orbit breakdown from existing results --- *)
detR1 = allResults[name]["details"]["R1"];
detR3 = allResults[name]["details"]["R3"];
Print["\n--- Sub-orbit breakdown RS ---"];
pairsR1 = Select[detR1, Length[#[[3]]] >= 2 &];
pairSub1 = {getSubtype[#[[3,1]]], getSubtype[#[[3,2]]], #[[4]]} & /@ pairsR1;
gr1 = GroupBy[pairSub1, {#[[1]], #[[2]]} &];
Do[
  entries = gr1[k];
  nsuc = Count[entries, _?(#[[3]] &)];
  ntot = Length[entries];
  Print["  ", k, "  n=", ntot, "  success: ", nsuc, "/", ntot]
, {k, Keys[gr1]}];

Print["\n--- Sub-orbit breakdown RZ ---"];
pairsR3 = Select[detR3, Length[#[[3]]] >= 2 &];
pairSub3 = {getSubtype[#[[3,1]]], getSubtype[#[[3,2]]], #[[4]]} & /@ pairsR3;
gr3 = GroupBy[pairSub3, {#[[1]], #[[2]]} &];
Do[
  entries = gr3[k];
  nsuc = Count[entries, _?(#[[3]] &)];
  ntot = Length[entries];
  Print["  ", k, "  n=", ntot, "  success: ", nsuc, "/", ntot]
, {k, Keys[gr3]}];

(* --- Face-type equivariance check --- *)
Print["\n--- Face-type equivariance check ---"];
rules = {"R1", "R3"};
Do[
  details = allResults[name]["details"][rule];
  validDetails = Select[details, Length[#[[3]]] >= 2 &];
  gr = GroupBy[validDetails,
    {faceType[[ #[[3, 1]] ]], faceType[[ #[[3, 2]] ]]} &];
  ruleOK = True;
  KeyValueMap[Function[{key, entries},
    successes = Count[entries, _?(#[[4]] &)];
    total     = Length[entries];
    If[successes != 0 && successes != total,
      ruleOK = False;
      Print["  [MIXED] ", rule, "  ", key, ": ", successes, "/", total]
    ]
  ], gr];
  If[ruleOK, Print["  [OK]   ", rule]]
, {rule, rules}];

(* --- 各 top でどの面が選ばれているか --- *)
Print["\n--- First few ordering sequences (RS) ---"];
Do[
  r = detR1[[i]];
  If[Length[r[[3]]] >= 3,
    Print["  top=", r[[3,1]], "(", getSubtype[r[[3,1]]], ")  ",
          "f2=", r[[3,2]], "(", getSubtype[r[[3,2]]], ")  ",
          "len=", Length[r[[3]]], "  ok=", r[[4]]]
  ]
, {i, 1, Min[10, Length[detR1]]}];
