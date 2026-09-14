#!/usr/bin/env wolframscript
(* net_separation_margin.wls
   How close does a face-rotation BFS net come to overlapping?

   The theorem says no two cells share interior volume, and the volume test
   returns exactly zero every time.  That is not the same as saying the
   verdict is robust: a referee is entitled to ask whether floating point
   could have hidden a small overlap.  This measures the margin.

   For each pair of cells not joined by a tree edge we compute the exact
   separating distance by the separating-axis construction: for two convex
   polyhedra the direction realising the minimum distance is a face normal
   of one of them or a cross product of one edge from each, so the maximum
   gap over that finite set of axes IS the distance.  Everything is dot
   products, which is also what an exact certificate would evaluate.

   Cells are congruent, so pairs whose centroids are further apart than
   twice the circumradius are separated by at least that excess; only the
   near pairs need the axis computation.

   Two numbers matter.  Pairs that TOUCH -- sharing a face, edge or vertex
   -- have distance exactly 0 and are common; they are not overlaps.  What
   must be far from zero is the smallest POSITIVE distance: if nothing sits
   between "touching" and "clearly apart", no rounding error can move a
   pair across the verdict.
*)
base = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
tol  = 10.^-9;
cutoff = 0.15;  (* pairs further apart than this are not resolved exactly *)

composeAff[{A2_, t2_}, {A1_, t1_}] := {A2 . A1, A2 . t1 + t2};
applyAff[{A_, t_}, x_] := A . x + t;
idAff4 = {IdentityMatrix[4], {0., 0., 0., 0.}};
rot4[e1_, e2_, th_] := IdentityMatrix[4] +
  (Cos[th]-1)(Outer[Times,e1,e1]+Outer[Times,e2,e2]) +
  Sin[th](Outer[Times,e2,e1]-Outer[Times,e1,e2]);

(* hull data once per cell, not once per pair: the hull is by far the
   most expensive part and each cell appears in many pairs *)
tri3 = Subsets[Range[4], {3}]; tri2 = Subsets[Range[4], {2}];
hullData[A_] := If[Length[A] == 4,
  {Table[Module[{q = A[[t]]}, Cross[q[[2]]-q[[1]], q[[3]]-q[[1]]]], {t, tri3}],
   Table[A[[t[[2]]]] - A[[t[[1]]]], {t, tri2}]},
  Module[{h = ConvexHullMesh[A]},
  {Select[Cross[#[[2]]-#[[1]], #[[3]]-#[[1]]] & /@
      (A[[#]] & /@ MeshCells[h,2][[All,1]]), Norm[#] > 10.^-9 &],
   Select[(#[[2]]-#[[1]]) & /@ (A[[#]] & /@ MeshCells[h,1][[All,1]]),
      Norm[#] > 10.^-9 &]}]];

(* Any axis gives a lower bound on the distance, so the face normals
   alone -- eight of them, against the dozens of edge-edge cross products
   -- settle every pair that is comfortably apart.  Only pairs the cheap
   bound leaves below the cutoff need the full axis set, and they are a
   small minority.  Axes are left unnormalised, the gaps divided by their
   lengths at the end. *)
gapAlong[axes_, A_, B_] := Module[{nrm, keep, ax, PA, PB},
  nrm = Sqrt[Total[axes^2, {2}]];
  keep = Flatten@Position[nrm, x_ /; x > 10.^-9];
  ax = axes[[keep]]; nrm = nrm[[keep]];
  PA = ax . Transpose[A]; PB = ax . Transpose[B];
  Max[MapThread[Max[#3 - #2, #1 - #4] &,
      {Min /@ PA, Max /@ PA, Min /@ PB, Max /@ PB}] / nrm]];

crossAxes[eA_, eB_] := Transpose[{
  Flatten[Outer[Times, eA[[All,2]], eB[[All,3]]] - Outer[Times, eA[[All,3]], eB[[All,2]]]],
  Flatten[Outer[Times, eA[[All,3]], eB[[All,1]]] - Outer[Times, eA[[All,1]], eB[[All,3]]]],
  Flatten[Outer[Times, eA[[All,1]], eB[[All,2]]] - Outer[Times, eA[[All,2]], eB[[All,1]]]]}];

satDistance[A_, dA_, B_, dB_] := Module[{g},
  g = gapAlong[Join[dA[[1]], dB[[1]]], A, B];
  If[g > cutoff, g,
     Max[g, gapAlong[crossAxes[dA[[2]], dB[[2]]], A, B]]]];

analyse[label_, file_, roots_] := Module[
  {raw, vers0, faces, cells, nC, f2c, adj, cellVerts, out = {}, R},
  raw = Get[base <> file];
  vers0 = N@raw[[1]]; faces = raw[[4]]; cells = raw[[6]]; nC = Length[cells];
  f2c = Merge[Flatten[Table[# -> i & /@ cells[[i]], {i, nC}]], Identity];
  adj = Table[{}, {nC}];
  Do[If[Length[f2c[k]] == 2,
     adj[[f2c[k][[1]]]] = Append[adj[[f2c[k][[1]]]], f2c[k][[2]]];
     adj[[f2c[k][[2]]]] = Append[adj[[f2c[k][[2]]]], f2c[k][[1]]]], {k, Keys[f2c]}];
  adj = Union /@ adj;
  cellVerts = Table[Union @@ faces[[cells[[i]]]], {i, nC}];

  Do[Module[{loc, cc, wHat = {0.,0.,0.,1.}, vv, ang, Rm, tfs, par, q, p,
             P, ctr, dm, close, ds},
    loc = (# - Mean[vers0]) & /@ vers0;
    cc = Mean[loc[[cellVerts[[root]]]]];
    vv = Which[Abs[(cc/Norm[cc]).wHat - 1] < 10.^-5, loc,
               Abs[(cc/Norm[cc]).wHat + 1] < 10.^-5, -loc,
               True, ang = ArcCos[Clip[cc.wHat/Norm[cc], {-1.,1.}]];
                     Rm = RotationMatrix[ang, {cc/Norm[cc], wHat}]; Rm.# & /@ loc];
    tfs = Association[root -> idAff4]; par = Association[root -> 0]; q = {root};
    While[q =!= {}, p = First[q]; q = Rest[q];
      Do[If[!KeyExistsQ[tfs, nb],
         Module[{fI, fv, fvT, cf, rel, u1, u2, v, nP, nCc, cosT, perp, e1, e2, th},
           fI = First[Intersection[cells[[p]], cells[[nb]]]];
           fv = vv[[faces[[fI]]]]; fvT = applyAff[tfs[p], #] & /@ fv; cf = Mean[fvT];
           rel = (# - cf) & /@ fvT; u1 = Normalize[rel[[1]]];
           u2 = Normalize[rel[[2]] - (rel[[2]].u1) u1];
           v = applyAff[tfs[p], Mean[vv[[cellVerts[[p]]]]]] - cf;
           nP = Normalize[v - (v.u1) u1 - (v.u2) u2];
           v = applyAff[tfs[p], Mean[vv[[cellVerts[[nb]]]]]] - cf;
           nCc = Normalize[v - (v.u1) u1 - (v.u2) u2];
           cosT = nCc.(-nP); e1 = nCc; perp = -nP - cosT nCc;
           e2 = If[Norm[perp] > 10.^-8, Normalize[perp],
                   Normalize@First@NullSpace[{u1,u2,nCc}]];
           th = ArcCos[Clip[cosT, {-1,1}]];
           tfs[nb] = composeAff[{rot4[e1,e2,th], cf - rot4[e1,e2,th].cf}, tfs[p]];
           par[nb] = p; AppendTo[q, nb]]], {nb, adj[[p]]}]];
    P = Table[(applyAff[tfs[i], #] & /@ vv[[cellVerts[[i]]]])[[All, 1;;3]], {i, nC}];
    ctr = Mean /@ P;
    R = Max[Table[Max[Norm /@ ((# - ctr[[i]]) & /@ P[[i]])], {i, nC}]];
    dm = DistanceMatrix[ctr];
    (* Pattern-matching Position over a 600x600 matrix costs more than the
       geometry it selects; this does the same with array arithmetic. *)
    close = Module[{tre = ConstantArray[0, {nC, nC}], mask},
      Do[If[par[i] > 0, tre[[i, par[i]]] = tre[[par[i], i]] = 1], {i, nC}];
      mask = UpperTriangularize[UnitStep[2 R + cutoff - dm] (1 - tre), 1];
      SparseArray[mask]["NonzeroPositions"]];
    Module[{hd = hullData /@ P},
      ds = Table[satDistance[P[[pr[[1]]]], hd[[pr[[1]]]],
                             P[[pr[[2]]]], hd[[pr[[2]]]]], {pr, close}]];
    AppendTo[out, {Length[close], Count[ds, x_ /; x < tol],
                   Min[Select[ds, # >= tol &] /. {} -> {Infinity}],
                   Min[ds], 2 R}]],
   {root, roots}];
  Module[{np = Total[out[[All,1]]], nt = Total[out[[All,2]]],
          mp = Min[out[[All,3]]], mn = Min[out[[All,4]]], d = out[[1,5]]},
   Print[StringPadRight[label, 10], " roots ", StringPadLeft[ToString@Length[roots],4],
     "   near pairs ", StringPadLeft[ToString@np, 8],
     "   touching ", StringPadLeft[ToString@nt, 7],
     "   min distance ", NumberForm[mn, {8,6}],
     "   smallest POSITIVE gap ", NumberForm[mp, {8,6}],
     "   (cell diameter ", NumberForm[d, {5,3}], ")"]]];

Print["distance between cells not joined by a tree edge, tree = T_r"];
analyse["5-cell",   "_4DData/f5.m",   Range[5]];
analyse["8-cell",   "_4DData/f8.m",   Range[8]];
analyse["16-cell",  "_4DData/f16.m",  Range[16]];
analyse["24-cell",  "_4DData/f24.m",  Range[24]];
analyse["120-cell", "_4DData/f120.m", Range[120]];
analyse["600-cell", "_4DData/f600.m", Range[600]];
