(* curvature_battery.wls -- 2026-08-19
   Tests whether any curvature-based scalar predicts net validity.
   Part 1 computes, for all 50 generated species, the full list of edge
   defects, the total curvature Sum(defect * length), and the cell-graph
   radius and diameter; it writes curvature_stats.mx.
   Part 2 ranks 13 candidate statistics by how well they isolate the
   three problem species (x3o4o3x ALL INVALID, x5o3x3x and x5x3o3x
   MIXED) out of 50.

   RESULT: none separates.  The best, totCurv / nC^(2/3), puts them at
   ranks 6, 8, 10 of 50 -- but the two lowest-scoring species are the
   grand antiprism and the snub 24-cell, both ALL VALID.  Worse, the
   accumulated statistics destroy the one clean signal the raw defect
   had: under totCurv/nC^(2/3) the six regular polytopes land at ranks
   34, 43, 17, 24, 35, 3, scrambling the difficulty order that the raw
   defect reproduced exactly.  Curvature scalars are not the invariant.

   Edge lengths are uniform to 1e-15 in the generated data, so
   Sum(defect*length) = Sum(defect). *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
Get[baseDir <> "wythoff_gen_all.mx"]; Get[baseDir <> "uniform4D_special.m"];
defectsL[verts_, edges_, faces_, cells_] := Module[
  {v = N[verts, MachinePrecision], nC = Length[cells], cellsOf, cv},
  cv = Table[Union @@ faces[[cells[[c]]]], {c, nC}];
  cellsOf = Lookup[GroupBy[Flatten[Table[{p, c}, {c, nC}, {p, cv[[c]]}], 1], First -> Last],
     Range[Length[v]], {}];
  Table[Module[{p = e[[1]], q = e[[2]], d, m, cs, len},
    d = v[[q]] - v[[p]]; len = Norm[d]; d = d/len; m = (v[[p]] + v[[q]])/2;
    cs = Intersection[cellsOf[[p]], cellsOf[[q]]];
    {2 Pi - Total@Table[Module[{fs, us},
       fs = Select[cells[[c]], SubsetQ[faces[[#]], {p, q}] &];
       us = Table[Module[{u = Mean[v[[faces[[f]]]]] - m}, u = u - (u . d) d; u/Norm[u]], {f, fs}];
       ArcCos[Clip[us[[1]] . us[[2]], {-1., 1.}]]], {c, cs}], len}], {e, edges}]];
res = {};
all = SortBy[Join[w4All, sp4All], Length[#["data"][[2]]] &];
Do[Module[{d = rec["data"], dl, dg, len, nC, adj, g, ecc},
   dl = defectsL[d[[1]], d[[2]], d[[4]], d[[6]]];
   dg = dl[[All, 1]]/Degree; len = dl[[All, 2]];
   nC = Length[d[[6]]];
   adj = Module[{f2c},
     f2c = GroupBy[Flatten@MapIndexed[Function[{cf,k},{# -> First[k]} & /@ cf], d[[6]]], First -> Last];
     Union /@ Lookup[GroupBy[Join[#, Reverse /@ #] &@ Select[Values[f2c], Length[#] === 2 &],
         First -> Last], Range[nC], {}]];
   g = Graph[Range[nC], Flatten[Table[UndirectedEdge[i,j],{i,nC},{j,Select[adj[[i]],#>i&]}]]];
   ecc = VertexEccentricity[g, #] & /@ Range[nC];
   AppendTo[res, <|"sym" -> rec["sym"], "nC" -> nC, "nE" -> Length[dg],
      "lenSpread" -> Max[len] - Min[len],
      "minD" -> Min[dg], "meanD" -> Mean[dg], "maxD" -> Max[dg],
      "totCurv" -> Total[dg len], "radius" -> Min[ecc], "diam" -> Max[ecc]|>];
   Print[rec["sym"], " done  minD=", NumberForm[Min[dg],{6,3}],
         "  totCurv=", NumberForm[Total[dg len],{8,1}],
         "  radius=", Min[ecc], "  diam=", Max[ecc],
         "  edge-length spread=", NumberForm[Max[len]-Min[len], {6,2}]]],
 {rec, all}];
DumpSave[baseDir <> "curvature_stats.mx", res];
Print["ALLDONE"];
baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
Get[baseDir <> "curvature_stats.mx"];
bad = {"x3o4o3x", "x5o3x3x", "x5x3o3x"};
stats = <|
  "minD"                -> (#minD &),
  "meanD"               -> (#meanD &),
  "totCurv"             -> (#totCurv &),
  "totCurv/nC"          -> (#totCurv/#nC &),
  "diam"                -> (#diam &),
  "radius"              -> (#radius &),
  "minD*diam"           -> (#minD #diam &),
  "minD*radius"         -> (#minD #radius &),
  "meanD*diam"          -> (#meanD #diam &),
  "meanD*radius"        -> (#meanD #radius &),
  "(totCurv/nC)*diam"   -> (#totCurv/#nC #diam &),
  "minD*nC^(1/3)"       -> (#minD #nC^(1/3) &),
  "totCurv/nC^(2/3)"    -> (#totCurv/#nC^(2/3) &)
|>;
n = Length[res];
Print["species = ", n, "   targets = ", bad, "\n"];
Print[StringPadRight["statistic", 20], StringPadLeft["ranks of the 3 (ascending)", 30],
      "   best-case separation"];
Do[Module[{f = stats[k], vals, ord, rk, mx},
   vals = f /@ res;
   ord = Ordering[vals];
   rk = Sort@Table[First@FirstPosition[ord, First@FirstPosition[res, r_ /; r["sym"] === b]],
        {b, bad}];
   mx = Max[rk]; mn = Min[rk];
   Print[StringPadRight[k, 20], StringPadLeft[ToString[rk], 30],
     "   bottom-", mx, " contains all 3", If[mx <= 6, "   <<< CANDIDATE", ""],
     "   |  top-", n - mn + 1]], {k, Keys[stats]}];
Print["\nA perfect separator would put the three at ranks {1,2,3} (or {48,49,50})."];
