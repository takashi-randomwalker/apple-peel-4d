(* ================================================================
   run_archimedean_faceup.m
   ================================================================
   アルキメデス13種を Face-up あり条件（p3lRunAll）で再計算し，
   R1 / R2 / R3 の成功数を集計する．
   結果は archimedean_faceup_results.mx に保存する．
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "peeling3DLoxo.m"];

(* --- アルキメデス13種 --- *)
(* 注意: R1 は peeling3DLoxo.m の p3lRunAll[vers, faces, "maxphi"] を使用する．
   以前ここにあったインライン版 p3lRunAllR1 は単一候補優先が未適用だったため削除した． *)
solidNames = {
  "TruncatedTetrahedron",
  "Cuboctahedron",
  "TruncatedCube",
  "TruncatedOctahedron",
  "Rhombicuboctahedron",
  "TruncatedCuboctahedron",
  "SnubCube",
  "Icosidodecahedron",
  "TruncatedDodecahedron",
  "TruncatedIcosahedron",
  "Rhombicosidodecahedron",
  "TruncatedIcosidodecahedron",
  "SnubDodecahedron"
};

rules = {"maxphi", "loxo", "maxz"};
ruleLabels = <|"maxphi" -> "R1", "loxo" -> "R2", "maxz" -> "R3"|>;

allResults = <||>;

Print[""];
Print[StringRepeat["=", 60]];
Print["  Archimedean Solids: Face-up Re-computation"];
Print[StringRepeat["=", 60]];

Do[
  {vers, faces} = p3lExtractPolyhedron[name];
  nf  = Length[faces];
  adj = p3lBuildAdj[faces];
  totalPairs = Total[Length /@ adj];

  Print["\n-- ", name, "  (F=", nf, ", pairs=", totalPairs, ") --"];

  solidRes = <||>;

  (* R1 *)
  t0 = AbsoluteTime[];
  r1 = p3lRunAll[vers, faces, "maxphi"];
  dt = AbsoluteTime[] - t0;
  n1 = Length[Select[r1, #[[4]]&]];
  Print["  R1: ", n1, "/", totalPairs,
        "  (", NumberForm[N[100 n1/totalPairs], {4,1}], "%)",
        "  ", NumberForm[dt, {4,2}], " s"];
  solidRes["R1"] = r1;

  (* R2 *)
  t0 = AbsoluteTime[];
  r2 = p3lRunAll[vers, faces, "loxo"];
  dt = AbsoluteTime[] - t0;
  n2 = Length[Select[r2, #[[4]]&]];
  Print["  R2: ", n2, "/", totalPairs,
        "  (", NumberForm[N[100 n2/totalPairs], {4,1}], "%)",
        "  ", NumberForm[dt, {4,2}], " s"];
  solidRes["R2"] = r2;

  (* R3 *)
  t0 = AbsoluteTime[];
  r3 = p3lRunAll[vers, faces, "maxz"];
  dt = AbsoluteTime[] - t0;
  n3 = Length[Select[r3, #[[4]]&]];
  Print["  R3: ", n3, "/", totalPairs,
        "  (", NumberForm[N[100 n3/totalPairs], {4,1}], "%)",
        "  ", NumberForm[dt, {4,2}], " s"];
  solidRes["R3"] = r3;

  allResults[name] = <|
    "nFaces" -> nf, "totalPairs" -> totalPairs,
    "R1" -> n1, "R2" -> n2, "R3" -> n3,
    "details" -> solidRes
  |>;,
  {name, solidNames}
];

(* --- 結果の保存 --- *)
outFile = scriptDir <> "archimedean_faceup_results.mx";
DumpSave[outFile, allResults];
Print["\nSaved: ", outFile];

(* --- サマリ表 --- *)
Print[""];
Print[StringRepeat["=", 60]];
Print["  Summary (Face-up condition)"];
Print[StringRepeat["=", 60]];
Print[StringPadLeft["Solid", 30], "  ",
      StringPadLeft["F", 4], "  ",
      StringPadLeft["Pairs", 6], "  ",
      StringPadLeft["R1", 6], "  ",
      StringPadLeft["R2", 6], "  ",
      StringPadLeft["R3", 6]];
Print[StringRepeat["-", 60]];
Do[
  r = allResults[name];
  Print[StringPadLeft[name, 30], "  ",
        StringPadLeft[ToString[r["nFaces"]], 4], "  ",
        StringPadLeft[ToString[r["totalPairs"]], 6], "  ",
        StringPadLeft[ToString[r["R1"]], 6], "  ",
        StringPadLeft[ToString[r["R2"]], 6], "  ",
        StringPadLeft[ToString[r["R3"]], 6]],
  {name, solidNames}
];

Print["\nDone."];
