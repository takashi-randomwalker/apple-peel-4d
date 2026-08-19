(* edge_angular_defect.wls -- 2026-08-19
   For every edge e of a 4-polytope,

       defect(e) = 2 Pi - sum over cells c containing e of
                          (dihedral angle of c at e).

   defect = 0 means the cells fit exactly around e, i.e. locally a
   face-to-face tiling of R^3.  Convexity forces defect > 0, so the
   defect measures how far the polytope is from a Euclidean honeycomb:
   {4,3,3} = 90 deg, {4,3,4} = 0 (cubic honeycomb), {4,3,5} < 0
   (hyperbolic).  It is also the holonomy of the unfolding: going once
   around e, the composed face-rotations return a rotation by defect(e),
   so a net is the developing map of a cone structure with cone angle
   2 Pi - defect along the edges.

   Regular six (single value each, all edges equivalent):
     5-cell 148.41, 8-cell 90.00, 16-cell 77.88,
     24-cell 31.59, 120-cell 10.30, 600-cell 7.36 degrees.
   That is exactly the order of how hard they are.

   BUT it does NOT predict validity for the uniform polytopes:
     x5x3x3x has min defect 0.736 and is ALL VALID, the same value as
     the MIXED x5x3o3x; and x3o4o3x is the hardest species yet sits at
     10.529, above many ALL VALID ones.  Falsified in both directions. *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
Get[baseDir <> "wythoff_gen_all.mx"]; Get[baseDir <> "uniform4D_special.m"];
defects[verts_, edges_, faces_, cells_] := Module[
  {v = N[verts, MachinePrecision], nC = Length[cells], cellsOf, cv},
  cv = Table[Union @@ faces[[cells[[c]]]], {c, nC}];
  cellsOf = Lookup[GroupBy[Flatten[Table[{p, c}, {c, nC}, {p, cv[[c]]}], 1], First -> Last],
     Range[Length[v]], {}];
  Table[Module[{p = e[[1]], q = e[[2]], d, m, cs},
    d = v[[q]] - v[[p]]; d = d/Norm[d]; m = (v[[p]] + v[[q]])/2;
    cs = Intersection[cellsOf[[p]], cellsOf[[q]]];
    2 Pi - Total@Table[Module[{fs, us},
       fs = Select[cells[[c]], SubsetQ[faces[[#]], {p, q}] &];
       us = Table[Module[{u = Mean[v[[faces[[f]]]]] - m}, u = u - (u . d) d; u/Norm[u]], {f, fs}];
       ArcCos[Clip[us[[1]] . us[[2]], {-1., 1.}]]], {c, cs}]], {e, edges}]];
all = SortBy[Join[w4All, sp4All], Length[#["data"][[2]]] &];
Do[Module[{d = rec["data"], dg},
   dg = defects[d[[1]], d[[2]], d[[4]], d[[6]]]/Degree;
   Print[StringPadRight[rec["sym"], 10], " C=", StringPadLeft[ToString[Length[d[[6]]]], 5],
     " E=", StringPadLeft[ToString[Length[d[[2]]]], 6],
     "   MIN defect ", NumberForm[Min[dg], {7, 3}],
     "   max ", NumberForm[Max[dg], {7, 3}]]], {rec, all}];
Print["ALLDONE"];
