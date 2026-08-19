(* check_tiebreak_orbits.wls -- 2026-08-19
   For each regular convex 4-polytope, root r=1: does every cell with more
   than one candidate BFS parent have those candidates in a single orbit of
   Stab(r,v)?  If yes, no equivariant tie-break can separate them.

   Result: yes for the 5-, 8-, 16- and 24-cell; NO for the 120-cell
   (50 of 71 cells split) and the 600-cell (84 of 183 split).

   NOTE: Graph[Range[n], edges] with an explicit vertex list is required.
   GraphAutomorphismGroup permutes VertexList indices, and Graph[edges]
   alone does not order its vertices as Range[n]. *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
Do[
 raw = Get[baseDir <> "_4DData/" <> f]; cellsP = raw[[6]]; n = Length[cellsP];
 adj = Table[Select[Range[n], Function[j, j =!= i &&
      Length[Intersection[cellsP[[i]], cellsP[[j]]]] >= 1]], {i, n}];
 edges = Flatten[Table[UndirectedEdge[i, j], {i, n}, {j, Select[adj[[i]], # > i &]}]];
 g = Graph[Range[n], edges];                    (* explicit vertex list *)
 Print[StringPadRight[f, 7], "  VertexList == Range[n] : ", VertexList[g] === Range[n],
       If[f === "f16.m", "   isomorphic to Q4 : " <> ToString[IsomorphicGraphQ[g, HypercubeGraph[4]]], ""]];
 els = GroupElements[GraphAutomorphismGroup[g]];
 (* sanity: every element must really be an automorphism *)
 okAut = AllTrue[els, Function[p, Sort[Sort /@ (List @@@ edges /. x_Integer :> PermutationReplace[x, p])] ===
                                  Sort[Sort /@ (List @@@ edges)]]];
 Print["         |Aut| = ", Length[els], "   all elements verified as automorphisms : ", okAut];
 r = 1; lay = <|r -> 0|>; q = {r};
 While[q =!= {}, v = First[q]; q = Rest[q];
  Do[If[! KeyExistsQ[lay, u], lay[u] = lay[v] + 1; AppendTo[q, u]], {u, adj[[v]]}]];
 nMulti = 0; nTrans = 0; bad = {};
 Do[cand = Select[adj[[v]], lay[#] == lay[v] - 1 &];
    If[Length[cand] > 1, nMulti++;
      fix = Select[els, PermutationReplace[r, #] == r && PermutationReplace[v, #] == v &];
      orb = FixedPoint[Union[#, Flatten[Table[PermutationReplace[x, p], {x, #}, {p, fix}]]] &,
                       {First[cand]}];
      If[SubsetQ[orb, cand], nTrans++, AppendTo[bad, {v, cand, Intersection[orb, cand]}]]],
  {v, Complement[Range[n], {r}]}];
 Print["         cells with >1 candidate parent : ", nMulti,
       "   candidates form ONE Stab(r,v)-orbit : ", nTrans,
       If[bad === {}, "   -> ALWAYS", "   -> failures: " <> ToString[bad]]],
 {f, {"f5.m", "f8.m", "f16.m", "f24.m", "f120.m", "f600.m"}}];
