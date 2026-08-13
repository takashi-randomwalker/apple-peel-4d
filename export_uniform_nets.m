(* =====================================================================
   export_uniform_nets.m
   Export the face-rotation BFS nets of the 64 convex uniform
   4-polytopes to uniform_nets/ , one net per polytope at root 1, plus
   three extra roots covering the exceptional cases.

   Two formats per net:
     stl/  binary STL - each cell written as its own closed triangulated
           shell.  For viewing and printing; the polygonal face
           structure is lost to triangulation.
     off/  OFF - vertices deduplicated by position, faces kept as
           polygons.  The cell decomposition, which OFF has no notion
           of, is preserved in "# cell" comment lines that standard
           viewers ignore.

   Also writes manifest.csv and README.md.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
outDir  = baseDir <> "uniform_nets/";
SetDirectory[baseDir];
Quiet[CreateDirectory /@ {outDir, outDir <> "stl", outDir <> "off"}];

Get[baseDir <> "wythoff_gen_all.mx"];
Get[baseDir <> "uniform4D_prisms.m"];   (* also loads uniform4D_special.m *)
Module[{src, cut},
  src = Import[baseDir <> "face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];
Get[baseDir <> "bfs_sat.m"];

(* ---------------------------------------------------------------- *)
(* Schlafli t-notation from the ring pattern                          *)
(* ---------------------------------------------------------------- *)
tNotation[sym_] := Module[{nodes, marks, rings},
  nodes = StringCases[sym, "x" | "o"];
  marks = StringCases[sym, DigitCharacter];
  rings = Flatten@Position[nodes, "x"] - 1;
  "t_{" <> StringRiffle[ToString /@ rings, ","] <> "}{" <>
   StringRiffle[marks, ","] <> "}"];

(* ---------------------------------------------------------------- *)
(* Net at a given root                                               *)
(* ---------------------------------------------------------------- *)
setup[rec_] := Module[{data, faceToCells, pairs},
  data   = rec["data"];
  versP  = N[data[[1]], MachinePrecision];
  facesP = data[[4]];
  cellsP = data[[6]];
  nCP    = Length[cellsP];
  versCenteredP = (# - Mean[versP]) & /@ versP;
  faceToCells = GroupBy[
    Flatten@MapIndexed[Function[{cf, k}, {# -> First[k]} & /@ cf], cellsP],
    First -> Last];
  pairs = Select[Values[faceToCells], Length[#] === 2 &];
  adjListP = Union /@ Lookup[
    GroupBy[Join[pairs, Reverse /@ pairs], First -> Last], Range[nCP], {}];
  adjSetP = Association@Flatten[
    Table[Sort[{i, j}] -> True, {i, nCP}, {j, adjListP[[i]]}], 1];
  satPrepCells[]];

netAt[r_] := Module[{tfsR},
  versP = satFaceUp[r, versCenteredP];
  tfsR  = satBfsUnfold[r];
  Table[Module[{tf = tfsR[i]},
     Take[#, 3] & /@ (versP[[satCellIdx[[i]]]] . Transpose[tf[[1]]] +
                      Threaded[tf[[2]]])], {i, nCP}]];

(* polygons of cell i, as lists of 3D points *)
cellFacePts[unf_, i_] := Module[{loc},
  loc = AssociationThread[satCellIdx[[i]] -> Range[Length[satCellIdx[[i]]]]];
  unf[[i]][[Lookup[loc, facesP[[#]]]]] & /@ cellsP[[i]]];

(* ---------------------------------------------------------------- *)
(* Binary STL                                                        *)
(* ---------------------------------------------------------------- *)
writeSTL[file_, unf_] := Module[{tris, str, n},
  tris = Join @@ Table[
    Join @@ (Function[poly,
       Table[poly[[{1, k, k + 1}]], {k, 2, Length[poly] - 1}]] /@
      cellFacePts[unf, i]), {i, nCP}];
  str = OpenWrite[file, BinaryFormat -> True];
  BinaryWrite[str, ConstantArray[0, 80], "Byte"];
  BinaryWrite[str, Length[tris], "UnsignedInteger32"];
  Do[
   n = Cross[t[[2]] - t[[1]], t[[3]] - t[[1]]];
   n = If[Norm[n] > 10^-12, n/Norm[n], {0., 0., 1.}];
   BinaryWrite[str, N@Join[n, t[[1]], t[[2]], t[[3]]], "Real32"];
   BinaryWrite[str, 0, "UnsignedInteger16"],
   {t, tris}];
  Close[str];
  Length[tris]];

(* ---------------------------------------------------------------- *)
(* OFF, polygons preserved, cell membership in comments              *)
(* ---------------------------------------------------------------- *)
writeOFF[file_, unf_] := Module[
  {grid = 10^-9, key, vmap, vbag, nv, verts, faces, cellOf, lines, idx},
  key[p_] := Round[p/grid];
  vmap = <||>; vbag = Internal`Bag[]; nv = 0;
  faces = Internal`Bag[]; cellOf = Internal`Bag[];
  Do[
   Do[
    idx = Table[
      With[{k = key[p]},
       If[! KeyExistsQ[vmap, k],
        Internal`StuffBag[vbag, p]; nv++; vmap[k] = nv];
       vmap[k]], {p, poly}];
    Internal`StuffBag[faces, idx]; Internal`StuffBag[cellOf, i],
    {poly, cellFacePts[unf, i]}],
   {i, nCP}];
  verts  = Internal`BagPart[vbag, All];
  faces  = Internal`BagPart[faces, All];
  cellOf = Internal`BagPart[cellOf, All];
  lines = Join[
    {"OFF",
     "# face-rotation BFS net; polygonal faces are preserved",
     "# the cell decomposition is listed below as comments",
     ToString[Length[verts]] <> " " <> ToString[Length[faces]] <> " 0"},
    (StringRiffle[ToString[CForm[N[#]]] & /@ #, " "]) & /@ verts,
    (StringRiffle[
       Prepend[ToString /@ (# - 1), ToString[Length[#]]], " "]) & /@ faces,
    Module[{byCell = GroupBy[Range[Length[cellOf]], cellOf[[#]] &]},
     Table["# cell " <> ToString[i] <> ": " <>
       StringRiffle[ToString /@ (Lookup[byCell, i, {}] - 1), " "],
      {i, nCP}]]];
  Export[file, StringRiffle[lines, "\n"], "Text"];
  {Length[verts], Length[faces]}];

(* ---------------------------------------------------------------- *)
(* What to export                                                    *)
(* ---------------------------------------------------------------- *)
w4Distinct = DeleteDuplicatesBy[w4All, {#["counts"], #["hist"]} &];

jobs = Join[
  Table[<|"rec" -> r, "sym" -> r["sym"], "name" -> tNotation[r["sym"]],
          "family" -> "Wythoffian", "root" -> 1|>, {r, w4Distinct}],
  Table[<|"rec" -> r, "sym" -> r["sym"], "name" -> r["name"],
          "family" -> "diminished 600-cell", "root" -> 1|>, {r, sp4All}],
  Table[<|"rec" -> r, "sym" -> StringReplace[r["name"], " " -> "_"],
          "name" -> r["name"], "family" -> "prism", "root" -> 1|>,
        {r, pr4All}],
  (* the exceptional cases, extra roots *)
  {<|"rec" -> SelectFirst[w4All, #["sym"] === "x3o4o3x" &],
     "sym" -> "x3o4o3x", "name" -> tNotation["x3o4o3x"],
     "family" -> "Wythoffian", "root" -> 47|>,
   <|"rec" -> SelectFirst[w4All, #["sym"] === "x5o3x3x" &],
     "sym" -> "x5o3x3x", "name" -> tNotation["x5o3x3x"],
     "family" -> "Wythoffian", "root" -> 0|>,
   <|"rec" -> SelectFirst[w4All, #["sym"] === "x5x3o3x" &],
     "sym" -> "x5x3o3x", "name" -> tNotation["x5x3o3x"],
     "family" -> "Wythoffian", "root" -> 748|>}];

(* root 0 means "the first failing root", resolved below *)
Get[baseDir <> "h4_bfs_allroots.mx"];
jobs = Replace[jobs,
  j_ /; j["root"] === 0 :> Module[{h},
    h = SelectFirst[h4All, #["sym"] === j["sym"] &];
    Append[j, "root" -> First@Select[h["results"], ! #[[4]] &][[1]]]],
  {1}];

Print["exporting ", Length[jobs], " nets to ", outDir];
Print[StringPadRight["file", 34], StringPadLeft["cells", 6],
      StringPadLeft["root", 6], "  verdict   ",
      StringPadLeft["tris", 8], StringPadLeft["STL MB", 9]];

rows = Table[
  Module[{j = jobs[[k]], base, chk, nt, off, sz},
   setup[j["rec"]];
   base = j["sym"] <> "_root" <> ToString[j["root"]];
   chk  = checkRootSAT[j["root"]];
   Module[{unf = netAt[j["root"]]},
    nt  = writeSTL[outDir <> "stl/" <> base <> ".stl", unf];
    off = writeOFF[outDir <> "off/" <> base <> ".off", unf]];
   sz = FileByteCount[outDir <> "stl/" <> base <> ".stl"]/1024.^2;
   Print[StringPadRight[base, 34], StringPadLeft[ToString[nCP], 6],
     StringPadLeft[ToString[j["root"]], 6], "  ",
     StringPadRight[If[chk[[4]], "valid", "SELF-INT"], 10],
     StringPadLeft[ToString[nt], 8], StringPadLeft[ToString[Round[sz, 0.01]], 9]];
   {base, j["sym"], j["name"], j["family"], nCP, j["root"],
    If[chk[[4]], "valid", "self-intersecting"], chk[[3]], off[[1]], off[[2]], nt}],
  {k, Length[jobs]}];

Export[outDir <> "manifest.csv",
  Prepend[rows, {"file", "coxeter_symbol", "name", "family", "cells",
                 "root", "verdict", "intersecting_pairs",
                 "net_vertices", "net_faces", "stl_triangles"}]];

Print["\ntotal STL ", Round[Total[FileByteCount /@
   FileNames["*.stl", outDir <> "stl"]]/1024.^2, 0.1], " MB"];
Print["total OFF ", Round[Total[FileByteCount /@
   FileNames["*.off", outDir <> "off"]]/1024.^2, 0.1], " MB"];
Print["manifest.csv written with ", Length[rows], " rows"];
