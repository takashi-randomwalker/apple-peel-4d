(* ================================================================
   check_snubdodec_fix.m
   ================================================================
   修正済み peeling3DLoxo.m で SnubDodecahedron を再計算し，
   旧版（archimedean_faceup_results.mx）と比較する．
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

(* サブタイプ定義 *)
pentFaces = Select[Range[nf], faceType[[#]] == 5 &];
triFaces  = Select[Range[nf], faceType[[#]] == 3 &];
triAdjPentCount = Association @ Table[
  f -> Length[Intersection[adj[[f]], pentFaces]], {f, triFaces}];
snubTri   = Select[triFaces, triAdjPentCount[#] >= 1 &];
gyrateTri = Select[triFaces, triAdjPentCount[#] == 0 &];
getSubtype[f_] := Which[
  faceType[[f]] == 5, "pent", MemberQ[snubTri, f], "snub", True, "gyrate"];

(* --- 旧版結果から face 42 の ペアの挙動 --- *)
Print["=== Old results (from .mx, old algorithm) ==="];
detOld = allResults[name]["details"]["R1"];
face42pairs = Select[detOld, #[[3,1]] == 42 &];
Print["Face 42 as F1: ", Length[face42pairs], " pairs"];
Do[
  r = face42pairs[[i]];
  Print["  F2=", r[[3,2]], "(", getSubtype[r[[3,2]]], ")  len=", Length[r[[3]]], "  ok=", r[[4]]]
, {i, Length[face42pairs]}];

(* --- 修正済みアルゴリズムで再計算 --- *)
Print["\n=== Re-computing with fixed algorithm ==="];
r1new = p3lRunAll[vers, faces, "maxphi"];
r3new = p3lRunAll[vers, faces, "maxz"];
n1new = Length[Select[r1new, #[[4]]&]];
n3new = Length[Select[r3new, #[[4]]&]];

n1old = allResults[name]["R1"];
n3old = allResults[name]["R3"];

Print["RS: old=", n1old, "/", totalPairs, "  new=", n1new, "/", totalPairs];
Print["RZ: old=", n3old, "/", totalPairs, "  new=", n3new, "/", totalPairs];

(* face 42 のペアを新版で確認 *)
face42new = Select[r1new, #[[3,1]] == 42 &];
Print["\nFace 42 as F1 (new algorithm):"];
Do[
  r = face42new[[i]];
  Print["  F2=", r[[3,2]], "(", getSubtype[r[[3,2]]], ")  len=", Length[r[[3]]], "  ok=", r[[4]]]
, {i, Length[face42new]}];

(* サブ軌道の変化確認 *)
pairsOld = Select[detOld, Length[#[[3]]] >= 2 &];
pairsNew = Select[r1new, Length[#[[3]]] >= 2 &];
subOld = GroupBy[{getSubtype[#[[3,1]]], getSubtype[#[[3,2]]], #[[4]]} & /@ pairsOld, {#[[1]], #[[2]]} &];
subNew = GroupBy[{getSubtype[#[[3,1]]], getSubtype[#[[3,2]]], #[[4]]} & /@ pairsNew, {#[[1]], #[[2]]} &];

Print["\n--- Sub-orbit comparison (RS) ---"];
allKeys = Union[Keys[subOld], Keys[subNew]];
Do[
  nOld = If[KeyExistsQ[subOld, k], Count[subOld[k], _?(#[[3]]&)], 0];
  nNew = If[KeyExistsQ[subNew, k], Count[subNew[k], _?(#[[3]]&)], 0];
  tot  = If[KeyExistsQ[subOld, k], Length[subOld[k]], Length[subNew[k]]];
  marker = If[nOld != nNew, " ← CHANGED", ""];
  Print["  ", k, "  n=", tot, "  old: ", nOld, "/", tot, "  new: ", nNew, "/", tot, marker]
, {k, allKeys}];
