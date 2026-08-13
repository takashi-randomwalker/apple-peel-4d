(* =====================================================================
   uniform4D_prisms.m
   The 17 prismatic convex uniform 4-polytopes: the prisms over the 5
   Platonic and 13 Archimedean solids, minus the cube prism, which is
   the tesseract and already counted among the regular ones.

   A prism over a uniform polyhedron P is a uniform 4-polytope provided
   the lateral edge length equals the edge length of P, so P is rescaled
   to unit edge length and the height is 1.

   Everything is combinatorial - no hull computation is needed:
     vertices  2n            (two copies of V(P))
     cells     F+2           (two copies of P, plus one prism per face)
     2-faces   2F + E        (two copies of each face, plus a rectangle
                              per edge)
     edges     2E + n        (two copies of each edge, plus the verticals)
   Euler: 2n-(2E+n)+(2F+E)-(F+2) = (n-E+F) - 2 = 0, using Euler for P.

   Incidence is then handed to sp4Incidence from uniform4D_special.m,
   which derives the ridges from cell pairs; C <= 64 here so its O(C^2)
   scan is free.
   ===================================================================== *)

ClearAll[pr4Solids, pr4Solid3D, pr4Prism, pr4Report];

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
Get[baseDir <> "uniform4D_special.m"];   (* for sp4Incidence, sp4Tol *)

(* the cube prism is the tesseract, so Cube is omitted *)
pr4Solids = {
  "Tetrahedron", "Octahedron", "Dodecahedron", "Icosahedron",
  "TruncatedTetrahedron", "Cuboctahedron", "TruncatedCube",
  "TruncatedOctahedron", "Rhombicuboctahedron", "TruncatedCuboctahedron",
  "SnubCube", "Icosidodecahedron", "TruncatedDodecahedron",
  "TruncatedIcosahedron", "Rhombicosidodecahedron",
  "TruncatedIcosidodecahedron", "SnubDodecahedron"};

(* vertices rescaled to unit edge length, faces as cyclic vertex lists *)
pr4Solid3D[name_] := Module[{v, f, e, len},
  v = N[PolyhedronData[name, "VertexCoordinates"], 30];
  f = PolyhedronData[name, "FaceIndices"];
  e = DeleteDuplicates[
    Sort /@ Flatten[Partition[#, 2, 1, 1] & /@ f, 1]];
  len = Norm[v[[e[[1, 1]]]] - v[[e[[1, 2]]]]];
  {v/len, f, e}];

(* 4D prism: bottom copy is 1..n, top copy is n+1..2n *)
pr4Prism[name_] := Module[{v, f, e, n, verts, cellsV},
  {v, f, e} = pr4Solid3D[name];
  n = Length[v];
  verts = Join[
    Append[#, -1/2] & /@ v,
    Append[#, +1/2] & /@ v];
  cellsV = Join[
    {Range[n], Range[n + 1, 2 n]},                 (* the two copies of P *)
    Table[Join[fc, fc + n], {fc, f}]];             (* one prism per face *)
  {verts, Sort /@ cellsV, {2 n, 2 Length[e] + n, 2 Length[f] + Length[e],
                           Length[f] + 2}}];

pr4Report[name_] := Module[{t, pr, r, v, e, ff, c, flags},
  {t, pr} = AbsoluteTiming[pr4Prism[name]];
  r = sp4Incidence[pr[[1]], pr[[2]]];
  {v, e, ff, c} = r["counts"];
  flags = {};
  If[{v, e, ff, c} =!= pr[[3]],
     AppendTo[flags, "counts!=" <> ToString[pr[[3]]]]];
  If[v - e + ff - c =!= 0, AppendTo[flags, "chi!=0"]];
  If[r["badRidges"] > 0, AppendTo[flags, "duplicate ridge"]];
  If[Length[r["degTally"]] =!= 1, AppendTo[flags, "degree not uniform"]];
  Print[StringPadRight[name, 28],
   StringPadLeft[ToString[v], 6], StringPadLeft[ToString[e], 7],
   StringPadLeft[ToString[ff], 7], StringPadLeft[ToString[c], 6],
   "  ", If[flags === {}, "OK      ", StringRiffle[flags, ","]],
   "cells ", r["hist"]];
  Join[r, <|"sym" -> name <> " prism", "name" -> name <> " prism",
            "ok" -> (flags === {})|>]];

Print["\n", StringRepeat["=", 92]];
Print["The 17 prismatic convex uniform 4-polytopes"];
Print[StringRepeat["=", 92]];
Print[StringPadRight["prism over", 28], StringPadLeft["V", 6],
      StringPadLeft["E", 7], StringPadLeft["F", 7], StringPadLeft["C", 6],
      "  status"];

pr4All = pr4Report /@ pr4Solids;

Print["\n", Count[pr4All, KeyValuePattern["ok" -> True]], " / ",
      Length[pr4All], " pass all checks"];
