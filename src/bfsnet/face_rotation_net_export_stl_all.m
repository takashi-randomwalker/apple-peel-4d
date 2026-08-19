(* face_rotation_net_export_stl_all.m
   Face-rotation BFS net STL export for all 6 regular convex 4-polytopes.
   Root cell = 1 for every polytope.
   Output: face_rotation_net_{name}_root1.stl
   Based on face_rotation_net_viz_all.m (reuses the same BFS + hull pipeline).
*)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";

(* ── 4D affine utilities ──────────────────────────────────────────────────── *)
composeAff[{A2_, t2_}, {A1_, t1_}] := {A2 . A1, A2 . t1 + t2};
applyAff[{A_, t_}, x_]             := A . x + t;
idAff4 = {IdentityMatrix[4], {0., 0., 0., 0.}};

rot4DInPlane[e1_, e2_, theta_] :=
  IdentityMatrix[4] +
  (Cos[theta] - 1)(Outer[Times, e1, e1] + Outer[Times, e2, e2]) +
  Sin[theta](Outer[Times, e2, e1] - Outer[Times, e1, e2]);

affRotAroundPt[cf_, e1_, e2_, theta_] :=
  With[{R = rot4DInPlane[N @ e1, N @ e2, N @ theta]}, {R, cf - R . cf}];

(* ── Per-edge unfolding transform ────────────────────────────────────────── *)
unfoldEdge[pIdx_, cIdx_, pTF_, versL_, facesL_, cellsL_] :=
  Module[{fIdx, fv, fvT, cf, rel, u1, u2, v, nP, nC, cosT, perp, e1, e2, theta},
    fIdx = First[Intersection[cellsL[[pIdx]], cellsL[[cIdx]]]];
    fv   = N @ versL[[ facesL[[fIdx]] ]];
    fvT  = applyAff[pTF, #] & /@ fv;
    cf   = Mean[fvT];
    rel  = (# - cf) & /@ fvT;
    u1   = Normalize[rel[[1]]];
    u2   = Normalize[rel[[2]] - (rel[[2]] . u1) u1];
    v    = applyAff[pTF, N @ Mean[versL[[ Union @@ facesL[[ cellsL[[pIdx]] ]] ]]]] - cf;
    nP   = Normalize[v - (v . u1) u1 - (v . u2) u2];
    v    = applyAff[pTF, N @ Mean[versL[[ Union @@ facesL[[ cellsL[[cIdx]] ]] ]]]] - cf;
    nC   = Normalize[v - (v . u1) u1 - (v . u2) u2];
    cosT  = nC . (-nP);
    e1    = nC;
    perp  = -nP - cosT nC;
    e2    = If[Norm[perp] > 10^-8, Normalize[perp],
               Normalize @ First @ NullSpace[{u1, u2, nC}]];
    theta = ArcCos[Clip[cosT, {-1, 1}]];
    composeAff[affRotAroundPt[cf, e1, e2, theta], pTF]
  ];

(* ── Polytope specs: {label, data file, STL output} ─────────────────────── *)
specs = {
  {"5-cell",   "_4DData/f5.m",   "face_rotation_net_5cell_root1.stl"},
  {"8-cell",   "_4DData/f8.m",   "face_rotation_net_8cell_root1.stl"},
  {"16-cell",  "_4DData/f16.m",  "face_rotation_net_16cell_root1.stl"},
  {"24-cell",  "_4DData/f24.m",  "face_rotation_net_24cell_root1.stl"},
  {"120-cell", "_4DData/f120.m", "face_rotation_net_120cell_root1.stl"},
  {"600-cell", "_4DData/f600.m", "face_rotation_net_600cell_root1.stl"}
};

(* ── Main loop ───────────────────────────────────────────────────────────── *)
Do[
  Block[{label, dataFile, stlFile, raw, vers, faces, cells, nC,
         adjList, locV, ccIdx, cc1, wHat, angle, R,
         tfs, bfsLayer, bfsQ, par, nb,
         unf4D, unf3D, hulls, tH, maxLayer, netViz, tTotal},

    {label, dataFile, stlFile} = spec;
    Print["\n", StringRepeat["-", 50]];
    Print[">> ", label];
    tTotal = AbsoluteTime[];

    (* 1. Load data *)
    raw   = Get[baseDir <> dataFile];
    vers  = N @ raw[[1]];
    faces = raw[[4]];
    cells = raw[[6]];
    nC    = Length[cells];
    Print["   cells=", nC, "  faces=", Length[faces],
          "  verts=", Length[vers]];

    (* 2. Adjacency list *)
    adjList = Table[
      Select[Range[nC],
        Function[j, j =!= i && Length[Intersection[cells[[i]], cells[[j]]]] >= 1]],
      {i, nC}];

    (* 3. Face-up rotation (root = cell 1) *)
    wHat  = {0., 0., 0., 1.};
    locV  = (# - Mean[vers]) & /@ vers;
    ccIdx = Union @@ faces[[ cells[[1]] ]];
    cc1   = Mean[locV[[ ccIdx ]]];
    vers  = Which[
      Abs[(cc1 / Norm[cc1]) . wHat - 1] < 0.0001, locV,
      Abs[(cc1 / Norm[cc1]) . wHat + 1] < 0.0001, -locV,
      True,
        angle = ArcCos[Clip[cc1 . wHat / Norm[cc1], {-1., 1.}]];
        R     = RotationMatrix[angle, {cc1 / Norm[cc1], wHat}];
        R . # & /@ locV];

    (* 4. BFS unfolding + layer tracking *)
    tfs      = Association[1 -> idAff4];
    bfsLayer = Association[1 -> 0];
    bfsQ     = {1};
    While[bfsQ =!= {},
      par = First[bfsQ]; bfsQ = Rest[bfsQ];
      Do[nb = neigh;
         If[!KeyExistsQ[tfs, nb],
           tfs[nb]      = unfoldEdge[par, nb, tfs[par], vers, faces, cells];
           bfsLayer[nb] = bfsLayer[par] + 1;
           AppendTo[bfsQ, nb]],
         {neigh, adjList[[par]]}]];
    maxLayer = Max[Values[bfsLayer]];
    Print["   BFS layers: 0 to ", maxLayer];

    (* 5. Project to 3D *)
    unf4D = Table[
      applyAff[tfs[i], #] & /@ vers[[ Union @@ faces[[ cells[[i]] ]] ]],
      {i, nC}];
    unf3D = unf4D[[All, All, 1 ;; 3]];
    Print["   w-spread: ",
          Max[Flatten[unf4D[[All, All, 4]]]] - Min[Flatten[unf4D[[All, All, 4]]]]];

    (* 6. Build convex hull meshes *)
    {tH, hulls} = AbsoluteTiming[
      Table[ConvexHullMesh[unf3D[[i]]], {i, nC}]];
    Print["   Hulls built in ", Round[tH, 0.1], " s"];

    (* 7. Build Graphics3D with BFS-layer gray shading *)
    netViz = Graphics3D[
      Table[{Opacity[0.72], EdgeForm[{GrayLevel[0.3], Thin}],
             GrayLevel[0.3 + 0.6 bfsLayer[i] / maxLayer],
             GraphicsComplex[
               unf3D[[i]],
               Polygon @ MeshCells[hulls[[i]], 2][[All, 1]]]},
            {i, nC}],
      Boxed           -> False,
      Lighting        -> "Neutral",
      PlotLabel       -> Style[label <> " face-rotation net  (BFS root = cell 1)",
                               14, Bold, Black],
      ViewPoint       -> {2.4, -2.0, 1.8},
      SphericalRegion -> True,
      Background      -> Black];

    (* 8. Export STL *)
    Export[baseDir <> stlFile, netViz];
    Print["   Saved: ", stlFile,
          "  (", Round[(AbsoluteTime[] - tTotal), 0.1], " s total)"]
  ],
  {spec, specs}];

Print["\n", StringRepeat["=", 50]];
Print["All 6 STL files exported."];
