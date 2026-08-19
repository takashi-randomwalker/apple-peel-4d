(* ================================================================
   check_mirror_snubcube.m
   ================================================================
   SnubCube の元の多面体と鏡像で，成功サブ軌道が変わるかを調べる．
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "peeling3DLoxo.m"];

name = "SnubCube";
{vers, faces} = p3lExtractPolyhedron[name];
nf    = Length[faces];
adj   = p3lBuildAdj[faces];
totalPairs = Total[Length /@ adj];
faceType = Length /@ faces;

(* サブタイプ *)
squareFaces   = Select[Range[nf], faceType[[#]] == 4 &];
triFaces      = Select[Range[nf], faceType[[#]] == 3 &];
triAdjSqCount = Association @ Table[
  f -> Length[Intersection[adj[[f]], squareFaces]], {f, triFaces}];
snubTri   = Select[triFaces, triAdjSqCount[#] == 1 &];
gyrateTri = Select[triFaces, triAdjSqCount[#] == 0 &];
getSubtype[f_] := Which[
  faceType[[f]] == 4, "sq",
  MemberQ[snubTri, f], "snub",
  True, "gyrate"
];

(* 元の多面体と鏡像 *)
versMirror = vers * ConstantArray[{-1., 1., 1.}, Length[vers]];

r1orig = p3lRunAll[vers,       faces, "maxphi"];
r1mir  = p3lRunAll[versMirror, faces, "maxphi"];
r3orig = p3lRunAll[vers,       faces, "maxz"];
r3mir  = p3lRunAll[versMirror, faces, "maxz"];

subBreakdown[results_, label_] := Module[{pairs, grouped},
  pairs   = Select[results, Length[#[[3]]] >= 2 &];
  grouped = GroupBy[{getSubtype[#[[3,1]]], getSubtype[#[[3,2]]], #[[4]]} & /@ pairs,
              {#[[1]], #[[2]]} &];
  Print["\n--- Sub-orbit breakdown: ", label, " ---"];
  Do[
    entries = grouped[k];
    nsuc = Count[entries, _?(#[[3]] &)];
    ntot = Length[entries];
    Print["  ", k, "  n=", ntot, "  success: ", nsuc, "/", ntot]
  , {k, Keys[grouped]}]
];

subBreakdown[r1orig, "RS original"];
subBreakdown[r1mir,  "RS mirror"];
subBreakdown[r3orig, "RZ original"];
subBreakdown[r3mir,  "RZ mirror"];

(* 成功ペア（F1,F2）のセットを比較 *)
Print["\n--- Successful (F1,F2) pair sets comparison ---"];
Do[
  {rule, rOrig, rMir} = triple;
  succOrig = Sort[#[[3,{1,2}]] & /@ Select[rOrig, #[[4]]&]];
  succMir  = Sort[#[[3,{1,2}]] & /@ Select[rMir,  #[[4]]&]];
  same = succOrig === succMir;
  nBoth     = Length[Intersection[succOrig, succMir]];
  nOnlyOrig = Length[Complement[succOrig, succMir]];
  nOnlyMir  = Length[Complement[succMir, succOrig]];
  Print[rule, ":  same pairs? ", same,
        "  both=", nBoth, "  only-orig=", nOnlyOrig, "  only-mirror=", nOnlyMir]
, {triple, {
    {"RS", r1orig, r1mir},
    {"RZ", r3orig, r3mir}
  }}];
