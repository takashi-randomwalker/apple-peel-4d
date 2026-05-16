(* ================================================================
   update_snubcube_results.m
   ================================================================
   archimedean_faceup_results.mx の SnubCube エントリを
   修正済み peeling3DLoxo.m (antipodal-case fix) で再計算・更新する．
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "peeling3DLoxo.m"];
Get[scriptDir <> "archimedean_faceup_results.mx"];

name = "SnubCube";
{vers, faces} = p3lExtractPolyhedron[name];
nf    = Length[faces];
adj   = p3lBuildAdj[faces];
totalPairs = Total[Length /@ adj];

Print["Re-computing ", name, " ..."];
r1 = p3lRunAll[vers, faces, "maxphi"];
r2 = p3lRunAll[vers, faces, "loxo"];
r3 = p3lRunAll[vers, faces, "maxz"];
n1 = Length[Select[r1, #[[4]]&]];
n2 = Length[Select[r2, #[[4]]&]];
n3 = Length[Select[r3, #[[4]]&]];

Print["  R1 (RS): ", n1, "/", totalPairs];
Print["  R2:      ", n2, "/", totalPairs];
Print["  R3 (RZ): ", n3, "/", totalPairs];

allResults[name] = <|
  "nFaces" -> nf, "totalPairs" -> totalPairs,
  "R1" -> n1, "R2" -> n2, "R3" -> n3,
  "details" -> <|"R1" -> r1, "R2" -> r2, "R3" -> r3|>
|>;

outFile = scriptDir <> "archimedean_faceup_results.mx";
DumpSave[outFile, allResults];
Print["Saved: ", outFile];
Print["Done."];
