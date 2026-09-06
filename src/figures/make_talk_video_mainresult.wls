#!/usr/bin/env wolframscript
(* make_talk_video_mainresult.wls
   Turntables of the 120-cell and the 600-cell nets, for the two large
   panels of Figure 1 on slide 6.

   The two are rendered as separate frame sequences and laid back into
   the figure by build_talk_video.sh, so that the slide keeps the layout
   it has: four small panels above, these two below with a white gap.
   One film, two panels turning together -- one click, not two.

   The style is that of src/bfsnet/face_rotation_net_viz_all.m, which is
   where the panels of 260612AllFaceRotation.pdf came from: grey by BFS
   depth, translucent, thin edges, and the same view point, so frame 0
   sits where the printed figure sits.  SphericalRegion keeps the net
   from breathing as the camera goes round.

   Args: <out dir> [frames]   frames default 240, i.e. 8 s at 30 fps.
*)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";

outDir = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]],
            baseDir <> "attic/outputs/frames_main"];
nFrames = If[Length[$ScriptCommandLine] >= 3,
             ToExpression[$ScriptCommandLine[[3]]], 240];

(* Rendered square and oversized.  SphericalRegion fits the sphere
   around the net, so the net itself covers only the inscribed part of
   the frame and comes out visibly smaller than in the printed figure;
   build_talk_video.sh crops back to what is drawn and scales the result
   into the 562x560 panel.  Cropping is safe here precisely because
   SphericalRegion holds the size fixed as the camera goes round. *)
imgSize = 800;
bgCol   = GrayLevel[0.13];

composeAff[{A2_, t2_}, {A1_, t1_}] := {A2 . A1, A2 . t1 + t2};
applyAff[{A_, t_}, x_]             := A . x + t;
idAff4 = {IdentityMatrix[4], {0., 0., 0., 0.}};

rot4DInPlane[e1_, e2_, theta_] :=
  IdentityMatrix[4] +
  (Cos[theta] - 1)(Outer[Times, e1, e1] + Outer[Times, e2, e2]) +
  Sin[theta](Outer[Times, e2, e1] - Outer[Times, e1, e2]);

affRotAroundPt[cf_, e1_, e2_, theta_] :=
  With[{R = rot4DInPlane[N @ e1, N @ e2, N @ theta]}, {R, cf - R . cf}];

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

(* Face-rotation BFS net of one polytope, as ready-to-draw primitives. *)
netPrimitives[dataFile_] :=
  Module[{raw, vers, faces, cells, nC, f2c, adj, tfs, lay, q, p,
          unf3D, hulls, mx},
    raw   = Get[baseDir <> dataFile];
    vers  = N @ raw[[1]];  faces = raw[[4]];  cells = raw[[6]];
    nC    = Length[cells];

    (* face -> its two cells, then cell -> neighbours; the O(C^2) scan
       of the older scripts is far too slow at 600 cells. *)
    f2c = Merge[Flatten[Table[# -> i & /@ cells[[i]], {i, nC}]], Identity];
    adj = Table[{}, {nC}];
    Do[If[Length[f2c[k]] == 2,
       adj[[ f2c[k][[1]] ]] = Append[adj[[ f2c[k][[1]] ]], f2c[k][[2]]];
       adj[[ f2c[k][[2]] ]] = Append[adj[[ f2c[k][[2]] ]], f2c[k][[1]]]],
     {k, Keys[f2c]}];
    adj = Sort /@ adj;

    Module[{wHat = {0., 0., 0., 1.}, locV, cc1, angle, R},
      locV = (# - Mean[vers]) & /@ vers;
      cc1  = Mean[locV[[ Union @@ faces[[ cells[[1]] ]] ]]];
      vers = Which[
        Abs[(cc1 / Norm[cc1]) . wHat - 1] < 0.0001, locV,
        Abs[(cc1 / Norm[cc1]) . wHat + 1] < 0.0001, -locV,
        True,
          angle = ArcCos[Clip[cc1 . wHat / Norm[cc1], {-1., 1.}]];
          R     = RotationMatrix[angle, {cc1 / Norm[cc1], wHat}];
          R . # & /@ locV]];

    tfs = Association[1 -> idAff4];  lay = Association[1 -> 0];  q = {1};
    While[q =!= {},
      p = First[q]; q = Rest[q];
      Do[If[!KeyExistsQ[tfs, nb],
          tfs[nb] = unfoldEdge[p, nb, tfs[p], vers, faces, cells];
          lay[nb] = lay[p] + 1;
          AppendTo[q, nb]],
        {nb, adj[[p]]}]];
    mx = Max[Values[lay]];

    unf3D = Table[
      (applyAff[tfs[i], #] & /@ vers[[ Union @@ faces[[ cells[[i]] ]] ]])[[All, 1 ;; 3]],
      {i, nC}];
    hulls = Table[ConvexHullMesh[unf3D[[i]]], {i, nC}];
    Print["   ", dataFile, ": ", nC, " cells, layers 0..", mx];
    Table[{Opacity[0.72], EdgeForm[{GrayLevel[0.3], Thin}],
           GrayLevel[0.3 + 0.6 lay[i] / mx],
           GraphicsComplex[unf3D[[i]],
             Polygon @ MeshCells[hulls[[i]], 2][[All, 1]]]},
          {i, nC}]
  ];

(* The view point of the printed figure, swept through a full turn. *)
r0   = Norm[{2.4, -2.0}];
phi0 = ArcTan[2.4, -2.0];

render[prims_, dir_, k_] :=
  Export[FileNameJoin[{dir, StringTemplate["f``.png"][IntegerString[k, 10, 4]]}],
    Graphics3D[prims,
      Boxed           -> False,
      Lighting        -> "Neutral",
      ImageSize       -> {imgSize, imgSize},
      ViewPoint       -> {r0 Cos[phi0 + 2 Pi k / nFrames],
                          r0 Sin[phi0 + 2 Pi k / nFrames], 1.8},
      ViewVertical    -> {0, 0, 1},
      SphericalRegion -> True,
      Background      -> bgCol],
    ImageSize -> {imgSize, imgSize}];

Do[
  Module[{dir = FileNameJoin[{outDir, spec[[1]]}], prims, t},
   If[!DirectoryQ[dir], CreateDirectory[dir, CreateIntermediateDirectories -> True]];
   Print[">> ", spec[[1]]];
   {t, prims} = AbsoluteTiming[netPrimitives[spec[[2]]]];
   Print["   built in ", t, " s"];
   Do[render[prims, dir, k];
      If[Mod[k, 40] == 0, Print["   frame ", k, "/", nFrames]],
     {k, 0, nFrames - 1}];
   Print["   wrote ", nFrames, " frames to ", dir]],
  {spec, {{"c120", "_4DData/f120.m"}, {"c600", "_4DData/f600.m"}}}];
