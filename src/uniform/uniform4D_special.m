(* =====================================================================
   uniform4D_special.m
   The two convex uniform 4-polytopes that Wythoff's construction on a
   ringed Coxeter diagram does not reach:

     snub 24-cell     s3s4o3o   V=96,  E=432, F=480, C=144
                                (120 tetrahedra + 24 icosahedra)
     grand antiprism            V=100, E=500, F=720, C=320
                                (300 tetrahedra + 20 pentagonal antiprisms)

   Both are diminishings of the 600-cell, and both share one construction:

     remove a set R of vertices, keep the convex hull of the rest.  Cells
     of the hull are
       (a) the 600-cell tetrahedra that use no removed vertex, plus
       (b) one cell per removed vertex v, namely the hull of v's
           remaining neighbors -- v's 12 neighbors form an icosahedron,
           and the cap is what survives of it.

     snub 24-cell     R = an inscribed 24-cell (24 vertices, pairwise
                          non-adjacent, so each cap keeps all 12
                          neighbors -> icosahedron)
     grand antiprism  R = two orthogonal decagonal rings (20 vertices;
                          each ring vertex loses its 2 ring neighbors,
                          so each cap keeps 10 -> pentagonal antiprism)

   The 600-cell is built here in standard icosian coordinates rather than
   taken from the Wythoff generator, because the inscribed 24-cell is
   then simply the 24 unit/half-integer quaternions and needs no search.
   Its cells are the 4-cliques of the edge graph (each vertex link is an
   icosahedron with 20 faces, so 120*20/4 = 600 cells).

   Output records match uniform4D_wythoff.m: "data" is the 6-part
   _4DData/f*.m list, plus "cellsV", "counts", "hist".
   ===================================================================== *)

ClearAll[sp4Tol, sp4Verts600, sp4Cells600, sp4Adj600, sp4Rings,
         sp4Diminish, sp4Incidence, sp4Report];

sp4Tol = 10^-8;

(* ---------------------------------------------------------------- *)
(* The 600-cell in standard icosian coordinates                      *)
(* ---------------------------------------------------------------- *)

sp4Verts600 := sp4Verts600 = Module[{phi, v8, v16, v96, base, evens},
  phi = N[GoldenRatio, 30];
  v8  = Flatten[Table[
     Table[If[k === i, s, 0], {k, 4}], {i, 4}, {s, {1, -1}}], 1];
  v16 = Tuples[{1/2, -1/2}, 4] // N[#, 30] &;
  base  = N[{0, 1, 1/phi, phi}/2, 30];
  evens = Select[Permutations[Range[4]], Signature[#] === 1 &];
  v96 = DeleteDuplicatesBy[
    Flatten[Table[base[[p]] * s, {p, evens}, {s, Tuples[{1, -1}, 4]}], 1],
    Round[#/sp4Tol] &];
  Join[N[v8, 30], v16, v96]];

(* adjacency: 600-cell edges are the pairs at cos 36 = phi/2 *)
sp4Adj600 := sp4Adj600 = Module[{v, n, ip, cosE},
  v = sp4Verts600; n = Length[v];
  cosE = N[GoldenRatio/2, 30];
  ip = v . Transpose[v];
  Table[Select[Range[n], Abs[ip[[i, #]] - cosE] < 10^-6 &], {i, n}]];

(* cells = 4-cliques of the edge graph *)
sp4Cells600 := sp4Cells600 = Module[{g},
  g = Graph[Range[Length[sp4Verts600]],
       UndirectedEdge @@@ Union[Sort /@ Flatten[
         MapIndexed[Function[{nb, k}, {First[k], #} & /@ nb],
                    sp4Adj600], 1]]];
  Sort /@ FindClique[g, {4}, All]];

(* ---------------------------------------------------------------- *)
(* Vertex subsets to remove                                          *)
(* ---------------------------------------------------------------- *)

(* the inscribed 24-cell: the 8 unit and 16 half-integer quaternions,
   which are the first 24 vertices by construction *)
sp4Removed24 := Range[24];

(* all decagonal rings: 10 vertices spanning a common 2-plane *)
sp4Rings := sp4Rings = Module[{v, n, out, comp, d1, d2, ring},
  v = sp4Verts600; n = Length[v];
  out = Reap[
    Do[
     comp = NullSpace[{v[[i]], v[[j]]}, Tolerance -> sp4Tol];
     If[Length[comp] === 2,
      d1 = v . comp[[1]]; d2 = v . comp[[2]];
      ring = Select[Range[n],
        Abs[d1[[#]]] < 10^-6 && Abs[d2[[#]]] < 10^-6 &];
      If[Length[ring] === 10, Sow[ring]]],
     {i, n - 1}, {j, i + 1, n}]][[2]];
  If[out === {}, {}, DeleteDuplicates[First[out]]]];

(* two rings whose planes are orthogonal complements *)
sp4RemovedGA := sp4RemovedGA = Module[{rs, v, pair},
  rs = sp4Rings; v = sp4Verts600;
  pair = SelectFirst[Subsets[rs, {2}],
    Max[Abs[v[[#[[1]]]] . Transpose[v[[#[[2]]]]]]] < 10^-6 &];
  If[pair === Missing["NotFound"], $Failed, Union @@ pair]];

(* ---------------------------------------------------------------- *)
(* Diminishing and incidence                                         *)
(* ---------------------------------------------------------------- *)

sp4Diminish[removed_] := Module[{keep, idx, surv, caps},
  keep = Complement[Range[Length[sp4Verts600]], removed];
  idx  = AssociationThread[keep -> Range[Length[keep]]];
  surv = Select[sp4Cells600, ContainsNone[#, removed] &];
  caps = Table[Complement[sp4Adj600[[v]], removed], {v, removed}];
  {sp4Verts600[[keep]], Sort /@ (Lookup[idx, #] & /@ Join[surv, caps])}];

(* generic convex-4-polytope incidence from vertices + cell vertex sets.
   C <= 320 here, so the O(C^2) ridge scan is cheap. *)
sp4Incidence[verts_, cellsV_] := Module[
  {nv, nc, ridges, faces, faceCells, cyc, edges, edgeMap, cellsF,
   faceAdj, vertAdj, badRidges},
  nv = Length[verts]; nc = Length[cellsV];

  ridges = Reap[
    Do[Module[{sh = Intersection[cellsV[[a]], cellsV[[b]]]},
      If[Length[sh] >= 3 &&
         MatrixRank[(# - verts[[First[sh]]]) & /@ verts[[sh]],
                    Tolerance -> sp4Tol] === 2,
       Sow[{sh, a, b}]]],
     {a, nc - 1}, {b, a + 1, nc}]][[2]];
  ridges = If[ridges === {}, {}, First[ridges]];

  faces     = ridges[[All, 1]];
  faceCells = ridges[[All, {2, 3}]];
  badRidges = Length[faces] - Length[DeleteDuplicates[faces]];

  cyc[ix_] := Module[{pts, o, rel, u1, u2, ang},
    pts = verts[[ix]]; o = Mean[pts];
    rel = (# - o) & /@ pts;
    u1  = Normalize[First[rel]];
    u2  = Normalize@First@Select[(# - (# . u1) u1) & /@ rel,
                                 Norm[#] > sp4Tol &];
    ang = ArcTan[# . u1, # . u2] & /@ rel;
    ix[[Ordering[ang]]]];
  faces = cyc /@ faces;

  edgeMap = GroupBy[
    Flatten[MapIndexed[
      Function[{f, k}, {Sort[#], First[k]} & /@ Partition[f, 2, 1, 1]],
      faces], 1], First -> Last];
  edges = Keys[edgeMap];

  cellsF = Lookup[
    GroupBy[Flatten@MapIndexed[
       Function[{ab, k}, {# -> First[k]} & /@ ab], faceCells],
     First -> Last], Range[nc], {}];
  faceAdj = Table[
    DeleteCases[
     Union @@ Lookup[edgeMap, Sort /@ Partition[faces[[k]], 2, 1, 1]], k],
    {k, Length[faces]}];
  vertAdj = Lookup[
    GroupBy[Join[edges, Reverse /@ edges], First -> Last], Range[nv], {}];

  <|"data"   -> {verts, edges, vertAdj, faces, faceAdj, cellsF},
    "cellsV" -> cellsV,
    "counts" -> {nv, Length[edges], Length[faces], nc},
    "hist"   -> Sort@Tally[Length /@ cellsV],
    "degTally"  -> Sort@Tally[Length /@ vertAdj],
    "badRidges" -> badRidges|>];

(* ---------------------------------------------------------------- *)
(* Report                                                            *)
(* ---------------------------------------------------------------- *)

sp4Report[name_, sym_, removed_, expect_] := Module[
  {t, dv, r, v, e, f, c, flags},
  {t, dv} = AbsoluteTiming[sp4Diminish[removed]];
  r = sp4Incidence @@ dv;
  {v, e, f, c} = r["counts"];
  flags = {};
  If[{v, e, f, c} =!= expect,
     AppendTo[flags, "counts!=" <> ToString[expect]]];
  If[v - e + f - c =!= 0, AppendTo[flags, "chi!=0"]];
  If[r["badRidges"] > 0, AppendTo[flags, "duplicate ridge"]];
  If[Length[r["degTally"]] =!= 1, AppendTo[flags, "degree not uniform"]];
  Print[StringPadRight[name, 18], StringPadRight[sym, 10],
   StringPadLeft[ToString[v], 6], StringPadLeft[ToString[e], 7],
   StringPadLeft[ToString[f], 7], StringPadLeft[ToString[c], 6],
   "  ", If[flags === {}, "OK      ", StringRiffle[flags, ","]],
   "cells ", r["hist"], "   deg ", r["degTally"]];
  Join[r, <|"sym" -> sym, "name" -> name, "ok" -> (flags === {})|>]];

(* ---------------------------------------------------------------- *)
(* Run                                                               *)
(* ---------------------------------------------------------------- *)

Print["600-cell scaffold: ", Length[sp4Verts600], " vertices, ",
      Length[sp4Cells600], " cells (4-cliques), vertex degree ",
      Union[Length /@ sp4Adj600]];
Print["decagonal rings found: ", Length[sp4Rings]];
Print["grand antiprism removal set: ", Length[sp4RemovedGA], " vertices"];
Print[""];
Print[StringPadRight["polytope", 18], StringPadRight["symbol", 10],
      StringPadLeft["V", 6], StringPadLeft["E", 7], StringPadLeft["F", 7],
      StringPadLeft["C", 6], "  status"];

sp4All = {
  sp4Report["snub 24-cell", "s3s4o3o", sp4Removed24, {96, 432, 480, 144}],
  sp4Report["grand antiprism", "gap", sp4RemovedGA, {100, 500, 720, 320}]};

Print["\n", Count[sp4All, KeyValuePattern["ok" -> True]], " / ",
      Length[sp4All], " pass all checks"];
