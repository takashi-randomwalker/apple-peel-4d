(* ================================================================
   check_mirror.m
   ================================================================
   各アルキメデス多面体について，元の多面体と鏡像（x 座標を反転）で
   Apple-Peel 結果が一致するかを検証する．

   鏡像の作り方: vers_mirror = vers * {-1, 1, 1}
   （yz 平面で反射．行列式 det = -1 なので向きが反転）

   面の頂点インデックスリストはそのまま流用できる
   （アルゴリズムが使うのは面の重心のみで，頂点順序は無関係）．

   期待:
     アンフィキラル多面体（鏡像が合同）: 元 = 鏡像
     キラル多面体（SnubCube, SnubDodecahedron）: 元 ≠ 鏡像 の可能性
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "peeling3DLoxo.m"];

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

(* キラル多面体（鏡像が非合同）*)
chiralSolids = {"SnubCube", "SnubDodecahedron"};

Print[""];
Print[StringRepeat["=", 70]];
Print["  Mirror symmetry check for Archimedean solids"];
Print["  Rules: RS (maxphi) and RZ (maxz)"];
Print[StringRepeat["=", 70]];

results = {};

Do[
  {vers, faces} = p3lExtractPolyhedron[name];
  nf         = Length[faces];
  adj        = p3lBuildAdj[faces];
  totalPairs = Total[Length /@ adj];

  (* 鏡像頂点座標 *)
  versMirror = vers * ConstantArray[{-1., 1., 1.}, Length[vers]];

  (* 元の多面体 *)
  r1orig = p3lRunAll[vers,       faces, "maxphi"];
  r3orig = p3lRunAll[vers,       faces, "maxz"];
  (* 鏡像 *)
  r1mir  = p3lRunAll[versMirror, faces, "maxphi"];
  r3mir  = p3lRunAll[versMirror, faces, "maxz"];

  n1orig = Length[Select[r1orig, #[[4]]&]];
  n3orig = Length[Select[r3orig, #[[4]]&]];
  n1mir  = Length[Select[r1mir,  #[[4]]&]];
  n3mir  = Length[Select[r3mir,  #[[4]]&]];

  chiral = MemberQ[chiralSolids, name];
  diff1 = n1orig != n1mir;
  diff3 = n3orig != n3mir;

  tag = Which[
    diff1 || diff3, " *** DIFFER ***",
    chiral,         " (chiral, but same)",
    True,           ""
  ];

  Print[name, If[chiral, " [CHIRAL]", ""], tag];
  Print["  RS: orig=", n1orig, "/", totalPairs, "  mirror=", n1mir, "/", totalPairs];
  Print["  RZ: orig=", n3orig, "/", totalPairs, "  mirror=", n3mir, "/", totalPairs];

  AppendTo[results,
    <|"name" -> name, "chiral" -> chiral,
      "total" -> totalPairs,
      "RS_orig" -> n1orig, "RS_mir" -> n1mir,
      "RZ_orig" -> n3orig, "RZ_mir" -> n3mir|>]

, {name, solidNames}];

Print[""];
Print[StringRepeat["=", 70]];
differ = Select[results, (#["RS_orig"] != #["RS_mir"] || #["RZ_orig"] != #["RZ_mir"])&];
If[Length[differ] == 0,
  Print["RESULT: All solids give identical results for original and mirror image."],
  Print["RESULT: ", Length[differ], " solid(s) give different results for mirror image:"];
  Do[
    Print["  ", r["name"],
          "  RS: ", r["RS_orig"], " vs ", r["RS_mir"],
          "  RZ: ", r["RZ_orig"], " vs ", r["RZ_mir"]]
  , {r, differ}]
];
Print[StringRepeat["=", 70]];
