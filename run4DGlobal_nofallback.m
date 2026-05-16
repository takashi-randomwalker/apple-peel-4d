(* ================================================================
   run4DGlobal_nofallback.m
   ================================================================
   peeling4Df4.m の noFB オプションで全6胞体・全C1を計算し，
   RZ/RS × with/without fallback の成功数を比較する．
   フルデータは保存せずサマリのみ ans4DGlobal_nofb.mx に保存．
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];
Get[scriptDir <> "peeling4Df4.m"];
dataDir = scriptDir <> "_4DData/";

polytopeNames = <|5->"5-cell", 8->"8-cell", 16->"16-cell",
                  24->"24-cell", 120->"120-cell", 600->"600-cell"|>;

loadPoly[n_] := Module[{raw},
  raw = Get[dataDir <> "f" <> ToString[n] <> ".m"];
  <|"vers"->raw[[1]], "faces"->raw[[4]], "cells"->raw[[6]]|>
];

summarize[ans_] := Module[{allBools, tc, fc, uniq, cl},
  allBools = Flatten[Map[Last, ans, {2}]];
  tc   = Count[allBools, True];
  fc   = Count[allBools, False];
  uniq = Length[DeleteDuplicates[
    Select[Flatten[ans, 1], Last[#]&][[All, 2]]]];
  cl   = Which[fc == 0, "Perfect", tc == 0, "Impossible", True, "Possible"];
  <|"true"->tc, "false"->fc, "unique"->uniq, "class"->cl|>
];

runOne[vers_, faces_, cells_, rule_, nofb_] := Module[{ans, nc = Length[cells]},
  ans = Table[peeling4Df4[vers, faces, cells, top, rule, nofb], {top, nc}];
  summarize[ans]
];

Print[""];
Print[StringRepeat["=", 70]];
Print["  Apple-Peel 4D: fallback あり/なし 比較  (RZ and RS)"];
Print[StringRepeat["=", 70]];

results4DNoFB = <||>;

Do[
  poly  = loadPoly[n];
  vers  = N[poly["vers"]]; faces = poly["faces"]; cells = poly["cells"];
  nc    = Length[cells];
  Print["\n--- ", polytopeNames[n], " (", nc, " cells) ---"];

  t0 = AbsoluteTime[];
  sRZw = runOne[vers, faces, cells, "RZ", False];
  Print["  RZ(w): ", Round[AbsoluteTime[]-t0, 0.1], " s  |  ",
        sRZw["class"], "  (True: ", sRZw["true"], ", Unique: ", sRZw["unique"], ")"];

  t0 = AbsoluteTime[];
  sRZn = runOne[vers, faces, cells, "RZ", True];
  Print["  RZ(n): ", Round[AbsoluteTime[]-t0, 0.1], " s  |  ",
        sRZn["class"], "  (True: ", sRZn["true"], ", Unique: ", sRZn["unique"], ")"];

  t0 = AbsoluteTime[];
  sRSw = runOne[vers, faces, cells, "RS", False];
  Print["  RS(w): ", Round[AbsoluteTime[]-t0, 0.1], " s  |  ",
        sRSw["class"], "  (True: ", sRSw["true"], ", Unique: ", sRSw["unique"], ")"];

  t0 = AbsoluteTime[];
  sRSn = runOne[vers, faces, cells, "RS", True];
  Print["  RS(n): ", Round[AbsoluteTime[]-t0, 0.1], " s  |  ",
        sRSn["class"], "  (True: ", sRSn["true"], ", Unique: ", sRSn["unique"], ")"];

  results4DNoFB[n] = <|
    "name" -> polytopeNames[n], "n" -> n,
    "RZ_w" -> sRZw, "RZ_n" -> sRZn,
    "RS_w" -> sRSw, "RS_n" -> sRSn
  |>,

  {n, {5, 8, 16, 24, 120, 600}}
];

outFile = scriptDir <> "ans4DGlobal_nofb.mx";
DumpSave[outFile, results4DNoFB];
Print["\nSaved: ", outFile];

Print[""];
Print[StringRepeat["=", 70]];
Print["  Summary"];
Print[StringRepeat["=", 70]];
Print[StringPadLeft["Polytope",  10], "  ",
      StringPadLeft["RZ(w)",     12], "  ",
      StringPadLeft["RZ(n)",     12], "  ",
      StringPadLeft["RS(w)",     12], "  ",
      StringPadLeft["RS(n)",     12]];
Print[StringRepeat["-", 62]];
Do[
  r = results4DNoFB[n];
  fmt[s_] := StringPadLeft[s["class"] <> "/" <> ToString[s["true"]], 12];
  Print[StringPadLeft[r["name"], 10], "  ",
        fmt[r["RZ_w"]], "  ", fmt[r["RZ_n"]], "  ",
        fmt[r["RS_w"]], "  ", fmt[r["RS_n"]]],
  {n, {5, 8, 16, 24, 120, 600}}
];
Print["Done."];
