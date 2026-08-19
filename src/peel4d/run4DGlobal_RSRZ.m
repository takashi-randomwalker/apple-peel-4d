(* ================================================================
   run4DGlobal_RSRZ.m
   ================================================================
   peeling4Df4.m（RS / RZ 両規則）を全6胞体・全 C1 に適用し比較する．
   結果を変数 results4D_RS, results4D_RZ に格納してコンソールに報告．
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];
Get[scriptDir <> "peeling4Df4.m"];

dataDir = scriptDir <> "_4DData/";

polytopeNames = <|
  5   -> "5-cell",
  8   -> "8-cell",
  16  -> "16-cell",
  24  -> "24-cell",
  120 -> "120-cell",
  600 -> "600-cell"
|>;

loadPoly[n_] := Module[{raw},
  raw = Get[dataDir <> "f" <> ToString[n] <> ".m"];
  <|"vers" -> raw[[1]], "faces" -> raw[[4]], "cells" -> raw[[6]]|>
];

summarize[ans_, nc_] := Module[
  {allBools, trueCount, falseCount, uniqOrders, classif},
  allBools   = Flatten[Map[Last, ans, {2}]];
  trueCount  = Count[allBools, True];
  falseCount = Count[allBools, False];
  uniqOrders = DeleteDuplicates[
    Select[Flatten[ans, 1], Last[#]&][[All, 2]]];
  classif = Which[falseCount == 0, "Perfect",
                  trueCount  == 0, "Impossible",
                  True,            "Possible"];
  <|"trueCount"      -> trueCount,
    "falseCount"     -> falseCount,
    "uniqueCount"    -> Length[uniqOrders],
    "classification" -> classif|>
];

Print[""];
Print[StringRepeat["=", 70]];
Print["  Apple-Peel 4D: RS vs RZ (all C1, all C2)"];
Print[StringRepeat["=", 70]];

results4D_RS = <||>;
results4D_RZ = <||>;

Do[
  poly = loadPoly[n];
  vers = poly["vers"]; faces = poly["faces"]; cells = poly["cells"];
  nc   = Length[cells];
  Print["\n--- ", polytopeNames[n], " (", nc, " cells) ---"];

  (* RZ *)
  t0  = AbsoluteTime[];
  ansRZ = Table[peeling4Df4[N[vers], faces, cells, top, "RZ"], {top, nc}];
  dtRZ  = Round[AbsoluteTime[] - t0, 0.1];
  sRZ   = summarize[ansRZ, nc];
  Print["  RZ: ", dtRZ, " s  |  ", sRZ["classification"],
        "  (True: ", sRZ["trueCount"],
        ", Unique: ", sRZ["uniqueCount"], ")"];

  (* RS *)
  t0  = AbsoluteTime[];
  ansRS = Table[peeling4Df4[N[vers], faces, cells, top, "RS"], {top, nc}];
  dtRS  = Round[AbsoluteTime[] - t0, 0.1];
  sRS   = summarize[ansRS, nc];
  Print["  RS: ", dtRS, " s  |  ", sRS["classification"],
        "  (True: ", sRS["trueCount"],
        ", Unique: ", sRS["uniqueCount"], ")"];

  results4D_RZ[n] = Join[sRZ, <|"name" -> polytopeNames[n]|>];
  results4D_RS[n] = Join[sRS, <|"name" -> polytopeNames[n]|>];

, {n, {5, 8, 16, 24, 120, 600}}];

(* ── サマリ表 ────────────────────────────────────── *)
Print[""];
Print[StringRepeat["=", 75]];
Print["  Summary: RS vs RZ"];
Print[StringRepeat["=", 75]];
Print[StringPadLeft["Polytope", 10], "  ",
      StringPadLeft["RZ Class", 10], "  ",
      StringPadLeft["RZ True",  8], "  ",
      StringPadLeft["RS Class", 10], "  ",
      StringPadLeft["RS True",  8]];
Print[StringRepeat["-", 55]];
Do[
  rz = results4D_RZ[n];
  rs = results4D_RS[n];
  Print[StringPadLeft[polytopeNames[n], 10], "  ",
        StringPadLeft[rz["classification"], 10], "  ",
        StringPadLeft[ToString[rz["trueCount"]], 8], "  ",
        StringPadLeft[rs["classification"], 10], "  ",
        StringPadLeft[ToString[rs["trueCount"]], 8]],
  {n, {5, 8, 16, 24, 120, 600}}
];

Print["\nDone."];
