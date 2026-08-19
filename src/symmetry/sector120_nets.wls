(* sector120_nets.wls -- 2026-08-19
   The 120-cell at maximum symmetry.  H = D_5 of order 10 (the
   stabilizer of a pentagonal face of the root dodecahedron; there are
   6 such axes).  Enumerates all 2048 H-invariant BFS trees, checks
   each for self-intersection, and ranks them by radius of gyration.

   Needs sector120_setup.mx, written by the setup stage (Stab(r), the
   layer map, the candidate parents and H).

   Result: all 2048 are valid; gyration 4.7224 .. 4.7811 over 563
   distinct values; the minimum is attained by exactly 2 trees.
   The queue-order tree we ship has symmetry order 1, so this is the
   one polytope where the symmetric construction gains something. *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
Module[{src, cut},
  src = Import[baseDir <> "src/bfsnet/face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];
Get[baseDir <> "src/bfsnet/bfs_sat.m"]; Get[baseDir <> "data/sector120_setup.mx"]; satPrepCells[];
img[pm_, x_] := PermutationReplace[x, pm];
orbits = DeleteDuplicates[Table[Union[Table[img[h, v], {h, H}]], {v, nCP}]];
reps = First /@ orbits;
fixedCand[v_] := With[{st = Select[H, img[#, v] === v &]},
   Select[cand[[v]], Function[c, AllTrue[st, img[#, c] === c &]]]];
nonRoot = DeleteCases[reps, r]; choices = fixedCand /@ nonRoot;
buildPar[sel_] := Module[{pp = <||>},
  MapThread[Function[{v, c}, Do[pp[img[h, v]] = img[h, c], {h, H}]], {nonRoot, sel}]; pp];
unfoldGiven[pm_] := Module[{tfs = <|r -> idAff4|>, order},
  order = SortBy[Complement[Range[nCP], {r}], lay[#] &];
  Do[tfs[v] = unfoldTFP[pm[v], v, tfs[pm[v]]], {v, order}];
  Table[Module[{tf = tfs[i]},
     Take[#,3] & /@ (versP[[satCellIdx[[i]]]] . Transpose[tf[[1]]] + Threaded[tf[[2]]])],
   {i, nCP}]];
gyr[u_] := Module[{pts = Flatten[u, 1], c}, c = Mean[pts]; Sqrt[Mean[Norm[# - c]^2 & /@ pts]]];
versP = satFaceUp[r, versCenteredP];
sels = Tuples[choices];
g = Table[gyr[unfoldGiven[buildPar[s]]], {s, sels}];
Print["H-invariant trees : ", Length[sels], "   |H| = ", Length[H]];
Print["radius of gyration : min ", NumberForm[Min[g], {9,6}],
      "   max ", NumberForm[Max[g], {9,6}],
      "   distinct values ", Length[Union[Round[g, 10.^-9]]]];
Print["how many attain the minimum : ", Count[g, x_ /; x < Min[g] + 10.^-9]];
best = First@Ordering[g, 1];
Print["best tree index : ", best];
uB = unfoldGiven[buildPar[sels[[best]]]];
(* validity of the best, by SAT *)
check[u_] := Module[{bMin, bMax, cent, radius, cd, ne},
  bMin = (Min /@ Transpose[#]) & /@ u; bMax = (Max /@ Transpose[#]) & /@ u;
  cent = Mean /@ u;
  radius = 2 Max[MapThread[Function[{pts, c}, Max[Norm[# - c] & /@ pts]], {u, cent}]];
  cd = satBBoxPairs[cent, bMin, bMax, radius];
  If[Length[cd] === 0, Return[0]];
  ne = Association@Table[c -> satNormalsEdges[u[[c]], satFaceTri[[c]], satEdgePair[[c]]],
     {c, Union @@ cd}];
  Length@Select[cd, Function[q, satOverlapQ[u[[q[[1]]]], ne[q[[1]]][[1]], ne[q[[1]]][[2]],
              u[[q[[2]]]], ne[q[[2]]][[1]], ne[q[[2]]][[2]]]]]];
Print["best net: overlaps = ", check[uB],
      "   R_gyration = ", NumberForm[Min[g], {9,6}],
      "   hull volume = ", NumberForm[RegionMeasure[ConvexHullMesh[Flatten[uB,1]]], {9,4}]];
uW = unfoldGiven[buildPar[sels[[Last@Ordering[g]]]]];
Print["worst invariant net: R_gyration = ", NumberForm[Max[g], {9,6}],
      "   hull volume = ", NumberForm[RegionMeasure[ConvexHullMesh[Flatten[uW,1]]], {9,4}],
      "   overlaps = ", check[uW]];
DumpSave[baseDir <> "data/sector120_best.mx", {uB, sels, best, g, H, nonRoot, choices}];
