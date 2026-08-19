(* ================================================================
   run4DPeeling.m
   ================================================================
   Runs peeling4Df on all 6 regular 4D polytopes and classifies
   each as Perfect / Possible / Impossible.

   Usage (in Mathematica):
     SetDirectory["/path/to/260324Peeling4D"];
     Get["peeling4Df.m"];
     Get["run4DPeeling.m"];
   ================================================================ *)

(* Path to the data folder -- NotebookDirectory[] for notebook, DirectoryName[$InputFileName] for .m *)
dataDir = If[ValueQ[$InputFileName], DirectoryName[$InputFileName], NotebookDirectory[]] <> "_4DData/";

(* Polytope names -- single line to avoid multi-line bracket issues *)
polytopeNames = Association[5 -> "5-cell (pentachoron)", 8 -> "8-cell (tesseract)", 16 -> "16-cell (hexadecachoron)", 24 -> "24-cell (icositetrachoron)", 120 -> "120-cell (hecatonicosachoron)", 600 -> "600-cell (hexacosichoron)"];

loadPolytope[n_] := Module[{dd, file, raw, vers, edgs, neis, faces, fneis, cells}, dd = If[ValueQ[dataDir], dataDir, NotebookDirectory[] <> "_4DData/"]; file = dd <> "f" <> ToString[n] <> ".m"; Print["Loading: ", file]; raw = Get[file]; If[!ListQ[raw] || Length[raw] != 6, Print["ERROR: Get returned ", Head[raw], " for ", file]; Return[$Failed]]; {vers, edgs, neis, faces, fneis, cells} = raw; Association["vers" -> vers, "faces" -> faces, "cells" -> cells, "n" -> n, "name" -> polytopeNames[n]]];

(* Classify one polytope as Perfect / Possible / Impossible *)
classifyPolytope[data_] := Module[
  {vers, faces, cells, n, ans, allBools, nCells, trueCount, falseCount, classification, t0},
  vers   = data["vers"];
  faces  = data["faces"];
  cells  = data["cells"];
  n      = data["n"];
  nCells = Length[cells];
  Print["Processing ", data["name"], " (", nCells, " cells) ..."];
  t0 = AbsoluteTime[];
  ans = Table[peeling4Df[N[vers], faces, cells, top], {top, nCells}];
  allBools  = Flatten[Map[Last, ans, {2}]];
  trueCount  = Count[allBools, True];
  falseCount = Count[allBools, False];
  classification = Which[falseCount == 0, "Perfect", trueCount == 0, "Impossible", True, "Possible"];
  Print["  Done in ", Round[AbsoluteTime[] - t0, 0.1], " s  |  ", classification, "  (True: ", trueCount, ", False: ", falseCount, ")"];
  Association["name" -> data["name"], "n" -> n, "classification" -> classification, "trueCount" -> trueCount, "falseCount" -> falseCount, "results" -> ans]
];

(* Main: run all polytopes and print summary table.
   Remove 120 and 600 from nList first if computation is too slow. *)
nList = {5, 8, 16, 24, 120, 600};

Print["======================================================"];
Print["  Apple-Peel Unfolding Classification for 4D Polytopes"];
Print["======================================================"];

results4D = Association @@ Table[
  Module[{data, res}, data = loadPolytope[n]; res = classifyPolytope[data]; n -> res],
  {n, nList}
];

Print["------ Summary ------"];
Print[TableForm[
  Prepend[
    Table[With[{r = results4D[n]}, {r["name"], r["classification"], r["trueCount"], r["falseCount"]}], {n, nList}],
    {"Polytope", "Classification", "True (F1,F2)", "False (F1,F2)"}
  ]
]];

(* Uncomment to save results:
DumpSave[NotebookDirectory[] <> "ans4D.mx", results4D]; *)

results4D
