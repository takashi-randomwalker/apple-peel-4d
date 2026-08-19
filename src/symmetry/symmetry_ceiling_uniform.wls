(* symmetry_ceiling_uniform.wls -- 2026-08-19
   The symmetry ceiling for all 64 convex uniform 4-polytopes: for one
   root per cell orbit, the largest subgroup H of Stab(r) admitting an
   H-invariant BFS tree.  Writes ceiling64.mx.  Runs in about 25 min.

   64 species, 170 cell orbits.  Ceiling per orbit:
     1 x28, 2 x10, 3 x21, 4 x19, 5 x17, 6 x42, 8 x17, 10 x14, 12 x1, 24 x1.
   28 of the 170 roots admit no symmetric net at all.
   Three species admit none from any root: the 16-cell, the 600-cell,
   and the runcinated 5-cell x3o3o3x.
   Two attain the full stabilizer at some root: the 5-cell (24) and the
   tetrahedron prism (12, at its orbit of four side cells).
   Only 18 of the 64 have a ceiling independent of the root; the grand
   antiprism ranges from 10 down to 1.

   IMPORTANT: the cell graph can have more automorphisms than the
   polytope has isometries.  For o3x3x3o the graph has 3840 but the
   polytope only 240, and taking the graph at face value gives a
   ceiling of 48 where the truth is 6.  The stabilizer is therefore
   filtered to the permutations induced by an orthogonal map on the
   cell centroids before the search. *)

(* Symmetry ceiling for all 64 convex uniform 4-polytopes.
   One root per cell orbit.  MAX |H| = largest subgroup of Stab(r)
   admitting an H-invariant BFS tree.  Each maximal H is also checked
   to be induced by an isometry of R^4, so the ceiling is geometric and
   not merely combinatorial. *)
baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
Get[baseDir <> "data/wythoff_gen_all.mx"];
Get[baseDir <> "src/uniform/uniform4D_special.m"];
Get[baseDir <> "src/uniform/uniform4D_prisms.m"];

(* drop the three cross-family duplicates by combinatorial signature *)
sig[rec_] := {Length[#[[1]]], Length[#[[2]]], Length[#[[4]]], Length[#[[6]]],
   Sort@Tally[Length /@ #[[6]]]} &[rec["data"]];
w45 = DeleteDuplicatesBy[w4All, sig];
Print["Wythoff species: ", Length[w4All], " -> ", Length[w45], " after dedup"];
all = Join[w45, sp4All, pr4All];
Print["total species: ", Length[all], "\n"];

analyse[rec_] := Module[
  {d = rec["data"], verts, cells, n, adj, g, aut, orbs, cents, out = {}},
  cells = d[[6]]; n = Length[cells];
  verts = N[d[[1]], MachinePrecision]; verts = (# - Mean[verts]) & /@ verts;
  cents = Table[Mean[verts[[Union @@ d[[4]][[cells[[c]]]]]]], {c, n}];
  adj = Module[{f2c},
    f2c = GroupBy[Flatten[Table[{f, c}, {c, n}, {f, cells[[c]]}], 1], First -> Last];
    Union /@ Lookup[GroupBy[Join[#, Reverse /@ #] &@ Select[Values[f2c], Length[#] === 2 &],
        First -> Last], Range[n], {}]];
  g = Graph[Range[n], Flatten[Table[UndirectedEdge[i, j], {i, n},
       {j, Select[adj[[i]], # > i &]}]]];
  aut = GraphAutomorphismGroup[g];
  orbs = SortBy[GroupOrbits[aut, Range[n]], Length];
  grow[gens_] := Module[{S = DeleteDuplicates[Prepend[gens, Range[n]]], new},
    While[True, new = DeleteDuplicates[Join[S, Flatten[Table[a[[b]], {a, S}, {b, gens}], 1]]];
     If[Length[new] === Length[S], Break[]]; S = new]; S];
  isometryQ[h_] := Module[{C2 = cents[[h]], A},
    A = LeastSquares[cents, C2];
    Max[Abs[A . Transpose[A] - IdentityMatrix[4]]] < 10.^-8 &&
    Max[Abs[cents . A - C2]] < 10.^-8];
  Do[Module[{r = Min[orb], stabA, lay, q, v, u, cand, feas, S1, frontier, maximal,
             next, grown, K2, H},
     r = Min[orb];
     stabA = grow[PermutationList[#, n] & /@ GroupGenerators[GroupStabilizer[aut, {r}]]];
     (* the cell graph can have more automorphisms than the polytope has
        isometries (o3x3x3o: 3840 vs 240), so keep only the realized ones *)
     stabA = Select[stabA, isometryQ];
     lay = <|r -> 0|>; q = {r};
     While[q =!= {}, v = First[q]; q = Rest[q];
      Do[If[! KeyExistsQ[lay, u], lay[u] = lay[v] + 1; AppendTo[q, u]], {u, Sort@adj[[v]]}]];
     cand = Table[If[v === r, {}, Select[adj[[v]], lay[#] == lay[v] - 1 &]], {v, n}];
     feas[H2_] := AllTrue[Complement[Range[n], {r}], Function[w,
        With[{st = Select[H2, #[[w]] === w &]},
          AnyTrue[cand[[w]], Function[c, AllTrue[st, #[[c]] === c &]]]]]];
     S1 = Select[stabA, feas[grow[{#}]] &];
     frontier = DeleteDuplicates[grow[{#}] & /@ S1]; maximal = {};
     While[frontier =!= {}, next = {};
      Do[grown = False;
         Do[K2 = grow[Join[K, {h}]];
            If[Length[K2] > Length[K] && feas[K2], grown = True; AppendTo[next, Sort[K2]]],
          {h, S1}];
         If[! grown, AppendTo[maximal, K]], {K, frontier}];
      frontier = DeleteDuplicates[next]];
     maximal = DeleteDuplicates[maximal];
     H = If[maximal === {}, {Range[n]}, First@MaximalBy[maximal, Length]];
     AppendTo[out, {r, Length[orb], Length[stabA], Length[H],
        AllTrue[H, isometryQ], Sort[Length /@ maximal]}]], {orb, orbs}];
  {GroupOrder[aut], Length /@ orbs, out}];

results = {};
Do[Module[{t0 = AbsoluteTime[], rec = all[[k]], a},
   a = analyse[rec];
   AppendTo[results, <|"sym" -> rec["sym"], "aut" -> a[[1]], "orbits" -> a[[2]],
      "perRoot" -> a[[3]]|>];
   Print[StringPadRight[ToString[k] <> ". " <> rec["sym"], 30],
     " C=", StringPadLeft[ToString[Length[rec["data"][[6]]]], 5],
     " |Aut|=", StringPadLeft[ToString[a[[1]]], 6],
     "  orbits ", a[[2]],
     "  ceilings ", a[[3]][[All, 4]],
     "  |Stab| ", a[[3]][[All, 3]],
     "  geometric ", Union[a[[3]][[All, 5]]],
     "  (", Round[AbsoluteTime[] - t0, 0.1], "s)"]],
 {k, Length[all]}];
DumpSave[baseDir <> "data/ceiling64.mx", results];
Print["\nALLDONE"];
