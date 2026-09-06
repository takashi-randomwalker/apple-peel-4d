#!/usr/bin/env wolframscript
(* make_talk_video_8cell.wls
   A turntable of the 8-cell face-rotation BFS net, for slide 2 of the
   JCDCG^3 talk.

   The geometry and the style are lifted unchanged from
   src/bfsnet/face_rotation_net_viz_all.m, so that frame 1 of the video
   and face_rotation_net_8cell.png are the same picture: the PNG stays
   the poster of the \movie, and a viewer without media support shows
   exactly what the slide showed before.  That is why the invisible
   PlotLabel is kept here too -- it reserves the top strip, and dropping
   it would recentre the net and break the match.

   SphericalRegion -> True is what keeps the net from breathing as the
   camera goes round; without it the bounding box is refitted every
   frame and the net pulses.

   Frames land in a directory given on the command line (default: the
   scratch directory below); build_talk_video.sh turns them into mp4.
*)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";

outDir = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]],
            baseDir <> "attic/outputs/frames_8cell"];
If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];

nFrames = 300;          (* 10 s at 30 fps *)
imgSize = 600;          (* same as the still *)
bgCol   = GrayLevel[0.13];

(* ── 4D affine utilities (from face_rotation_net_viz_all.m) ─────────────── *)
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

(* ── Load the 8-cell and unfold from root cell 1 ────────────────────────── *)
raw   = Get[baseDir <> "_4DData/f8.m"];
vers  = N @ raw[[1]];
faces = raw[[4]];
cells = raw[[6]];
nCell = Length[cells];

adjList = Table[
  Select[Range[nCell],
    Function[j, j =!= i && Length[Intersection[cells[[i]], cells[[j]]]] >= 1]],
  {i, nCell}];

wHat  = {0., 0., 0., 1.};
locV  = (# - Mean[vers]) & /@ vers;
cc1   = Mean[locV[[ Union @@ faces[[ cells[[1]] ]] ]]];
vers  = Which[
  Abs[(cc1 / Norm[cc1]) . wHat - 1] < 0.0001, locV,
  Abs[(cc1 / Norm[cc1]) . wHat + 1] < 0.0001, -locV,
  True,
    Module[{angle, R},
      angle = ArcCos[Clip[cc1 . wHat / Norm[cc1], {-1., 1.}]];
      R     = RotationMatrix[angle, {cc1 / Norm[cc1], wHat}];
      R . # & /@ locV]];

tfs      = Association[1 -> idAff4];
bfsLayer = Association[1 -> 0];
bfsQ     = {1};
While[bfsQ =!= {},
  Module[{par = First[bfsQ]},
    bfsQ = Rest[bfsQ];
    Do[If[!KeyExistsQ[tfs, nb],
        tfs[nb]      = unfoldEdge[par, nb, tfs[par], vers, faces, cells];
        bfsLayer[nb] = bfsLayer[par] + 1;
        AppendTo[bfsQ, nb]],
      {nb, adjList[[par]]}]]];
maxLayer = Max[Values[bfsLayer]];

unf4D = Table[
  applyAff[tfs[i], #] & /@ vers[[ Union @@ faces[[ cells[[i]] ]] ]], {i, nCell}];
unf3D = unf4D[[All, All, 1 ;; 3]];
Print["cells=", nCell, "  layers 0..", maxLayer,
      "  w-spread=", Max[Flatten[unf4D[[All, All, 4]]]] -
                     Min[Flatten[unf4D[[All, All, 4]]]]];

hulls = Table[ConvexHullMesh[unf3D[[i]]], {i, nCell}];

prims = Table[{Opacity[0.72], EdgeForm[{GrayLevel[0.3], Thin}],
               GrayLevel[0.3 + 0.6 bfsLayer[i] / maxLayer],
               GraphicsComplex[unf3D[[i]],
                 Polygon @ MeshCells[hulls[[i]], 2][[All, 1]]]},
              {i, nCell}];

(* ── Turntable about the vertical axis ─────────────────────────────────── *)
(* The still uses ViewPoint {2.4, -2.0, 1.8}; keep its radius and height
   and sweep the azimuth through a full turn, so the video starts and
   ends on exactly the still's camera. *)
r0   = Norm[{2.4, -2.0}];
phi0 = ArcTan[2.4, -2.0];

frame[k_] := Graphics3D[prims,
  Boxed           -> False,
  Lighting        -> "Neutral",
  ImageSize       -> {imgSize, imgSize},
  PlotLabel       -> Style["8-cell face-rotation net  (BFS root = cell 1)",
                           14, Bold, bgCol],
  ViewPoint       -> {r0 Cos[phi0 + 2 Pi k / nFrames],
                      r0 Sin[phi0 + 2 Pi k / nFrames], 1.8},
  ViewVertical    -> {0, 0, 1},
  SphericalRegion -> True,
  Background      -> bgCol];

Do[
  Export[FileNameJoin[{outDir, StringTemplate["f``.png"][IntegerString[k, 10, 4]]}],
         frame[k], ImageSize -> {imgSize, imgSize}];
  If[Mod[k, 50] == 0, Print["  frame ", k, "/", nFrames]],
  {k, 0, nFrames - 1}];

Print["Wrote ", nFrames, " frames to ", outDir];
