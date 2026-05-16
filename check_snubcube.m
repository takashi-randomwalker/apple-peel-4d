(* ================================================================
   SnubCube 詳細解析 (変数名を安全な名前に修正)
   ================================================================ *)
scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "peeling3DLoxo.m"];
Get[scriptDir <> "archimedean_faceup_results.mx"];

{vers, faces} = p3lExtractPolyhedron["SnubCube"];
adj           = p3lBuildAdj[faces];
faceType      = Length /@ faces;
nf            = Length[faces];

Print["SnubCube: ", nf, " faces, face-type counts: ", Counts[faceType]];

squareFaces   = Select[Range[nf], faceType[[#]] == 4 &];
triFaces      = Select[Range[nf], faceType[[#]] == 3 &];

(* 各三角形が正方形に辺接する回数 *)
triAdjSqCount = Table[
  f -> Length[Intersection[adj[[f]], squareFaces]],
  {f, triFaces}
] // Association;

snubTri   = Select[triFaces, triAdjSqCount[#] == 1 &];  (* 正方形に1辺接 *)
gyrateTri = Select[triFaces, triAdjSqCount[#] == 0 &];  (* 正方形に0辺接 *)

Print["\nTriangle subtypes:  snub=", Length[snubTri],
      "  gyrate=", Length[gyrateTri]];

(* 各面のサブタイプ *)
getSubtype[f_] := Which[
  faceType[[f]] == 4, "sq",
  MemberQ[snubTri, f],    "snub",
  MemberQ[gyrateTri, f],  "gyrate"
];

(* --- RS 結果の解析 --- *)
detR1 = allResults["SnubCube"]["details"]["R1"];

(* order[[1]]=F1面インデックス，order[[2]]=F2面インデックス *)
pairs = Select[detR1, Length[#[[3]]] >= 2 &];
Print["\nTotal valid pairs: ", Length[pairs]];

(* (F1-subtype, F2-subtype) ごとに集計 *)
pairWithSubtype = {getSubtype[#[[3,1]]], getSubtype[#[[3,2]]], #[[4]]} & /@ pairs;
grouped = GroupBy[pairWithSubtype, {#[[1]], #[[2]]} &];
Print["\n--- Sub-orbit breakdown (RS) ---"];
Do[
  entries = grouped[k];
  nsuc    = Count[entries, _?(#[[3]] &)];
  ntot    = Length[entries];
  Print["  ", k, "  n=", ntot, "  success: ", nsuc, "/", ntot]
, {k, Keys[grouped]}];

(* RZ も同様 *)
detR3 = allResults["SnubCube"]["details"]["R3"];
pairsR3 = Select[detR3, Length[#[[3]]] >= 2 &];
pairWithSubtypeR3 = {getSubtype[#[[3,1]]], getSubtype[#[[3,2]]], #[[4]]} & /@ pairsR3;
groupedR3 = GroupBy[pairWithSubtypeR3, {#[[1]], #[[2]]} &];
Print["\n--- Sub-orbit breakdown (RZ) ---"];
Do[
  entries = groupedR3[k];
  nsuc    = Count[entries, _?(#[[3]] &)];
  ntot    = Length[entries];
  Print["  ", k, "  n=", ntot, "  success: ", nsuc, "/", ntot]
, {k, Keys[groupedR3]}];

(* --- (sq, snub) 失敗ペアの詳細 --- *)
Print["\n--- Failed (sq, snub-triangle) pairs ---"];
failPairs = Select[pairs,
  getSubtype[#[[3,1]]] == "sq" &&
  getSubtype[#[[3,2]]] == "snub" &&
  !#[[4]] &];
If[Length[failPairs] == 0,
  Print["  None (all sq-snub pairs succeed)"],
  Do[
    f1i = p[[3, 1]]; f2i = p[[3, 2]];
    Print["  F1=", f1i, "  F2=", f2i,
          "  c(F1)=", N[Mean[vers[[ faces[[f1i]] ]]], 3],
          "  c(F2)=", N[Mean[vers[[ faces[[f2i]] ]]], 3]]
  , {p, failPairs}]
];
