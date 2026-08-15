(* =====================================================================
   export_uniform_data.m
   Write the combinatorial data of the 64 convex uniform 4-polytopes to
   _4DData_uniform/ , one file per polytope, in a 5-element format:

     [[1]] verts       vertex coordinates in R^4, machine precision
     [[2]] edges       vertex-index pairs
     [[3]] faces       2-faces as vertex-index lists, in CYCLIC order
     [[4]] cellsByFace cells as face-index lists
     [[5]] cellAdj     for each cell, the cells sharing a 2-face

   This is _4DData/f*.m with two fields dropped and one added.  The
   dropped ones, vertex adjacency and face adjacency, are never read by
   face_rotation_net_all4D_v2.m and cost 8.6 MB across the 45 Wythoffian
   members; the added one, cell adjacency, is what the unfolding
   actually needs and costs 1.3 MB.  See _4DData_uniform/README.md.

   Coordinates are written at machine precision because the pipeline
   applies N[...] on load; the 30-digit working precision of the
   generator is not meaningful downstream.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
outDir  = baseDir <> "_4DData_uniform/";
SetDirectory[baseDir];
Quiet[CreateDirectory[outDir]];

Get[baseDir <> "wythoff_gen_all.mx"];
Get[baseDir <> "uniform4D_prisms.m"];   (* also loads uniform4D_special.m *)

tNotation[sym_] := Module[{nodes, marks, rings},
  nodes = StringCases[sym, "x" | "o"];
  marks = StringCases[sym, DigitCharacter];
  If[Length[nodes] =!= 4 || Length[marks] =!= 3, Return[sym]];
  rings = Flatten@Position[nodes, "x"] - 1;
  "t_{" <> StringRiffle[ToString /@ rings, ","] <> "}{" <>
   StringRiffle[marks, ","] <> "}"];

knownNames = <|
  "x3o3o3o" -> "5-cell", "o3o3o3x" -> "5-cell",
  "x4o3o3o" -> "tesseract", "o4o3o3x" -> "16-cell",
  "x3o4o3o" -> "24-cell", "o3o4o3x" -> "24-cell",
  "x5o3o3o" -> "120-cell", "o5o3o3x" -> "600-cell",
  "x3o4o3x" -> "runcinated 24-cell",
  "x5o3x3x" -> "runcitruncated 600-cell",
  "x5x3o3x" -> "runcitruncated 120-cell",
  "s3s4o3o" -> "snub 24-cell", "gap" -> "grand antiprism"|>;

cellAdjacency[cellsF_] := Module[{ftc, prs},
  ftc = GroupBy[
    Flatten@MapIndexed[Function[{cf, k}, {# -> First[k]} & /@ cf], cellsF],
    First -> Last];
  prs = Select[Values[ftc], Length[#] === 2 &];
  Union /@ Lookup[GroupBy[Join[prs, Reverse /@ prs], First -> Last],
                  Range[Length[cellsF]], {}]];

w4Distinct = DeleteDuplicatesBy[w4All, {#["counts"], #["hist"]} &];

jobs = Join[
  Table[<|"rec" -> r, "base" -> r["sym"], "sym" -> r["sym"],
          "family" -> "Wythoffian"|>, {r, w4Distinct}],
  Table[<|"rec" -> r, "base" -> r["sym"], "sym" -> r["sym"],
          "family" -> "diminished 600-cell"|>, {r, sp4All}],
  Table[<|"rec" -> r, "base" -> StringReplace[r["name"], " " -> "_"],
          "sym" -> StringReplace[r["name"], " " -> "_"],
          "family" -> "prism"|>, {r, pr4All}]];

Print["writing ", Length[jobs], " files to _4DData_uniform/"];
Print[StringPadRight["file", 32], StringPadLeft["V", 7],
      StringPadLeft["E", 7], StringPadLeft["F", 7], StringPadLeft["C", 6],
      StringPadLeft["MB", 8]];

rows = Table[
  Module[{j = jobs[[k]], d, cadj, out, file, sz, v, e, f, c, nm},
   d    = j["rec"]["data"];
   cadj = cellAdjacency[d[[6]]];
   out  = {N[d[[1]]], d[[2]], d[[4]], d[[6]], cadj};
   file = outDir <> j["base"] <> ".m";
   Put[out, file];
   sz = FileByteCount[file]/1024.^2;
   {v, e, f, c} = j["rec"]["counts"];
   nm = Lookup[knownNames, j["sym"], tNotation[j["sym"]]];
   Print[StringPadRight[j["base"] <> ".m", 32],
     StringPadLeft[ToString[v], 7], StringPadLeft[ToString[e], 7],
     StringPadLeft[ToString[f], 7], StringPadLeft[ToString[c], 6],
     StringPadLeft[ToString[Round[sz, 0.01]], 8]];
   {j["base"] <> ".m", j["sym"], nm, j["family"], v, e, f, c,
    Round[Mean[N[Length /@ cadj]], 0.01], Round[sz, 0.001]}],
  {k, Length[jobs]}];

Export[outDir <> "index.csv",
  Prepend[rows, {"file", "coxeter_symbol", "name", "family",
                 "vertices", "edges", "faces", "cells",
                 "mean_cell_degree", "size_MB"}]];

Print["\ntotal ", Round[Total[FileByteCount /@ FileNames["*.m", outDir]]/1024.^2, 0.1],
      " MB in ", Length[FileNames["*.m", outDir]], " files"];

(* --- verify a round trip on one large and one small case ---------- *)
Print["\nround-trip check:"];
Do[Module[{raw = Get[outDir <> f <> ".m"], v, e, fa, cf, ca, ok},
   {v, e, fa, cf, ca} = raw;
   ok = And[
     Length[v] === Length[DeleteDuplicates[Round[v/10^-9]]],
     (* Euler *)
     Length[v] - Length[e] + Length[fa] - Length[cf] === 0,
     (* cell adjacency is symmetric and matches the shared-face rule *)
     And @@ Table[And @@ (MemberQ[ca[[#]], i] & /@ ca[[i]]), {i, Length[cf]}],
     ca === Table[Select[Range[Length[cf]],
        # =!= i && Intersection[cf[[i]], cf[[#]]] =!= {} &], {i, Length[cf]}]];
   Print["  ", StringPadRight[f, 12], If[ok, "OK", "*** FAILED ***"],
         "   V-E+F-C=", Length[v] - Length[e] + Length[fa] - Length[cf],
         "  mean cell degree ", Round[Mean[N[Length /@ ca]], 0.01]]],
 {f, {"o3o3o3x", "x5x3x3x", "gap", "SnubDodecahedron_prism"}}];
