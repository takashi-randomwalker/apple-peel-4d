(* symmetry_ceiling.wls -- 2026-08-19
   For each regular convex 4-polytope with root r = 1, find the largest
   subgroup H <= Stab(r) that admits an H-invariant BFS spanning tree.
   By the equivariance of face-rotation unfolding, the resulting net is
   then H-symmetric, so this is the ceiling on net symmetry.

   H-invariant tree exists  <=>  every cell v has a candidate parent
   fixed by Stab_H(v).  Feasible subgroups are downward closed, so
   climbing the lattice from feasible cyclic subgroups is exact.

   Usage: wolframscript -file symmetry_ceiling.wls '{"f5.m","f8.m"}'

   Result:  5-cell 24/24, 8-cell 8/48, 16-cell 1/24, 24-cell 6/48,
            120-cell 10/120, 600-cell 1/24.
   For cube / octahedron / dodecahedron cells the ceiling equals the
   stabilizer of one 2-face of the root cell; tetrahedral cells (16-,
   600-cell) collapse to the trivial group. *)

(* maximal H <= Stab(r) admitting an H-invariant BFS tree, r = 1.
   The feasible subgroups are downward closed, so climbing the lattice
   from cyclic subgroups finds every maximal one exactly. *)
baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
img[p_, x_] := PermutationReplace[x, p];
files = ToExpression[$ScriptCommandLine[[2]]];

Do[
 t0 = AbsoluteTime[];
 raw = Get[baseDir <> "_4DData/" <> f]; cellsP = raw[[6]]; n = Length[cellsP];
 adj = Table[Select[Range[n], Function[j, j =!= i &&
      Length[Intersection[cellsP[[i]], cellsP[[j]]]] >= 1]], {i, n}];
 g = Graph[Range[n], Flatten[Table[UndirectedEdge[i, j], {i, n},
       {j, Select[adj[[i]], # > i &]}]]];
 grp = GraphAutomorphismGroup[g];
 r = 1;
 stabR = GroupElements[GroupStabilizer[grp, {r}]];
 lay = <|r -> 0|>; q = {r};
 While[q =!= {}, v = First[q]; q = Rest[q];
  Do[If[! KeyExistsQ[lay, u], lay[u] = lay[v] + 1; AppendTo[q, u]], {u, Sort@adj[[v]]}]];
 cand = Table[If[v === r, {}, Select[adj[[v]], lay[#] == lay[v] - 1 &]], {v, n}];

 close[gens_] := Sort@GroupElements[PermutationGroup[gens]];
 feasibleQ[H_] := AllTrue[Complement[Range[n], {r}], Function[v,
    With[{st = Select[H, img[#, v] === v &]},
      AnyTrue[cand[[v]], Function[c, AllTrue[st, img[#, c] === c &]]]]]];

 S1 = Select[stabR, feasibleQ[close[{#}]] &];
 (* climb: start from feasible cyclic subgroups, keep adding generators *)
 frontier = DeleteDuplicates[close[{#}] & /@ S1];
 maximal = {};
 While[frontier =!= {},
  next = {};
  Do[grown = False;
     Do[K2 = close[Join[K, {h}]];
        If[Length[K2] > Length[K] && feasibleQ[K2],
           grown = True; AppendTo[next, K2]], {h, S1}];
     If[! grown, AppendTo[maximal, K]], {K, frontier}];
  frontier = DeleteDuplicates[next]];
 maximal = DeleteDuplicates[maximal];
 best = If[maximal === {}, 1, Max[Length /@ maximal]];

 Print[StringPadRight[f, 8],
   "  cells=", StringPadLeft[ToString[n], 4],
   "  |Stab(r)|=", StringPadLeft[ToString[Length[stabR]], 5],
   "  individually-feasible elements=", StringPadLeft[ToString[Length[S1]], 5],
   "  MAX |H| = ", StringPadLeft[ToString[best], 4],
   "  maximal subgroup orders = ", Sort[Length /@ maximal],
   "  (", Round[AbsoluteTime[] - t0, 0.1], " s)"],
 {f, files}];
