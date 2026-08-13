(* face_rotation_net_all4D_v2.m
   Face-rotation BFS net check for ALL 6 regular 4-polytopes.
   v2: Uses strict interior intersection (RegionMeasure > eps)
       instead of RegionMember (which falsely counts boundary contact).

   Corrected method:
     Two cells "intersect" iff RegionMeasure[RegionIntersection[h1,h2]] > eps
   This ignores shared edges/vertices (boundary contact) and detects
   only true interior overlaps.
*)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";

(* ================================================================
   Polytope-independent 4D affine utilities (unchanged from v1)
   ================================================================ *)
composeAff[{A2_,t2_},{A1_,t1_}] := {A2.A1, A2.t1+t2};
applyAff[{A_,t_},x_]            := A.x+t;
idAff4 = {IdentityMatrix[4], {0.,0.,0.,0.}};

rot4DInPlane[e1_,e2_,theta_] :=
  IdentityMatrix[4] +
  (Cos[theta]-1)(Outer[Times,e1,e1]+Outer[Times,e2,e2]) +
  Sin[theta](Outer[Times,e2,e1]-Outer[Times,e1,e2]);

affRotAroundPt[cf_,e1_,e2_,theta_] :=
  With[{R=rot4DInPlane[N@e1,N@e2,N@theta]}, {R, cf-R.cf}];

(* ================================================================
   Generic data-dependent functions (use globals versP,facesP,cellsP)
   ================================================================ *)
cVertsIdxP[i_] := Union @@ facesP[[ cellsP[[i]] ]];
cVertsP[i_]    := versP[[ cVertsIdxP[i] ]];
cCentroidP[i_] := Mean @ N @ cVertsP[i];

sharedFaceVertsP[i_,j_] :=
  versP[[ facesP[[ First[Intersection[cellsP[[i]], cellsP[[j]]]] ]] ]];

adjQP[i_,j_] := KeyExistsQ[adjSetP, Sort[{i,j}]];

faceUpForRootP[r_, locV_] := Module[
  {ccIdx,cc,wHat={0.,0.,0.,1.},angle,R},
  ccIdx = Union@@facesP[[cellsP[[r]]]];
  cc    = Mean[locV[[ccIdx]]];
  Which[
    Abs[(cc/Norm[cc]).wHat-1]<0.0001, locV,
    Abs[(cc/Norm[cc]).wHat+1]<0.0001, -locV,
    True,
      angle=ArcCos[Clip[cc.wHat/Norm[cc],{-1.,1.}]];
      RotationMatrix[angle,{cc/Norm[cc],wHat}].#&/@locV]];

unfoldTFP[pIdx_,cIdx_,pTF_] := Module[
  {fv,fvT,cf,u1,u2,rel,nPar,nChd,v,perp,e1,e2,cosT,theta},
  fv  =N@sharedFaceVertsP[pIdx,cIdx];
  fvT =applyAff[pTF,#]&/@fv;
  cf  =Mean[fvT];
  rel =(#-cf)&/@fvT;
  u1  =Normalize[rel[[1]]];
  u2  =Normalize[rel[[2]]-(rel[[2]].u1)u1];
  v   =applyAff[pTF,N@cCentroidP[pIdx]]-cf;
  nPar=Normalize[v-(v.u1)u1-(v.u2)u2];
  v   =applyAff[pTF,N@cCentroidP[cIdx]]-cf;
  nChd=Normalize[v-(v.u1)u1-(v.u2)u2];
  cosT=nChd.(-nPar);
  e1  =nChd;
  perp=-nPar-cosT*nChd;
  e2  =If[Norm[perp]>10^-8,Normalize[perp],
          Normalize@First@NullSpace[{u1,u2,nChd}]];
  theta=ArcCos[Clip[cosT,{-1,1}]];
  composeAff[affRotAroundPt[cf,e1,e2,theta],pTF]];

bfsUnfoldP[root_] := Module[
  {q={root},vis={root},tfs=Association[root->idAff4],par,nb},
  While[q=!={},
    par=First[q];q=Rest[q];
    Do[nb=neigh;
       If[!MemberQ[vis,nb],
         AppendTo[vis,nb];AppendTo[q,nb];
         tfs[nb]=unfoldTFP[par,nb,tfs[par]]],
       {neigh,adjListP[[par]]}]];
  tfs];

(* ================================================================
   v2 intersection check: strict interior via RegionMeasure
   Two cells truly intersect iff their 3D hull intersection has
   positive volume (> eps).  Boundary contact (shared face/edge/vertex)
   gives RegionMeasure = 0 and is correctly ignored.
   ================================================================ *)
epsVol = 10^-8;  (* volume threshold: adjust if needed *)

trueIntersectQ[hull1_, hull2_] :=
  Quiet[RegionMeasure[RegionIntersection[hull1, hull2]] > epsVol,
        {RegionMeasure::noReg, RegionIntersection::reg}];

(* ================================================================
   Check one root (v2: strict interior only)
   ================================================================ *)
checkRootV2[r_] := Module[
  {tfsR,unf4D,unf3D,wSpread,bMinR,bMaxR,bboxList,nBBox,hulls,trueHits},

  versP   = faceUpForRootP[r, versCenteredP];
  tfsR    = bfsUnfoldP[r];
  unf4D   = Table[applyAff[tfsR[i],#]&/@N@cVertsP[i], {i,nCP}];
  wSpread = (Max[#]-Min[#])&@Flatten[unf4D[[All,All,4]]];
  unf3D   = unf4D[[All,All,1;;3]];

  bMinR    = (Min/@Transpose[#])&/@unf3D;
  bMaxR    = (Max/@Transpose[#])&/@unf3D;
  bboxList = Reap[
    Do[If[!adjQP[i,j] &&
          And@@Thread[bMinR[[i]]<=bMaxR[[j]]] &&
          And@@Thread[bMinR[[j]]<=bMaxR[[i]]],
         Sow[{i,j}]],
       {i,nCP},{j,i+1,nCP}]][[2]];
  bboxList = If[bboxList==={},{},First[bboxList]];
  nBBox    = Length[bboxList];

  If[nBBox==0, Return[{r,0,0,True}]];

  hulls    = Table[ConvexHullMesh[unf3D[[i]]], {i,nCP}];
  trueHits = Select[bboxList, Function[p,
    trueIntersectQ[hulls[[p[[1]]]], hulls[[p[[2]]]]]]];

  {r, nBBox, Length[trueHits], Length[trueHits]==0}
];

(* ================================================================
   Run one polytope
   ================================================================ *)
runPolytopeV2[name_, dataFile_, dualNote_] := Module[
  {raw, tLoad, tAdj, tAll, results, nValid},

  Print["\n", StringRepeat["-",55]];
  Print["Polytope: ", name, "  (", dualNote, ")"];
  Print[StringRepeat["-",55]];

  {tLoad, raw} = AbsoluteTiming[Get[baseDir<>"_4DData/"<>dataFile]];
  versP  = raw[[1]]; facesP = raw[[4]]; cellsP = raw[[6]];
  nCP    = Length[cellsP];
  versCenteredP = (#-Mean[N@versP])&/@N@versP;
  Print["  vertices=",Length[versP],"  faces=",Length[facesP],
        "  cells=",nCP,"  (loaded in ",tLoad," s)"];

  {tAdj, adjListP} = AbsoluteTiming[
    Table[Select[Range[nCP],
      Function[j, j=!=i &&
        Length[Intersection[cellsP[[i]],cellsP[[j]]]]>=1]],
      {i,nCP}]];
  adjSetP = Association@Flatten[
    Table[Sort[{i,j}]->True,{i,nCP},{j,adjListP[[i]]}],1];
  Print["  adj degree: min=",Min[Length/@adjListP],
        "  max=",Max[Length/@adjListP],
        "  (built in ",tAdj," s)"];

  {tAll, results} = AbsoluteTiming[
    Table[checkRootV2[r], {r,1,nCP}]];
  nValid = Length[Select[results, #[[4]]&]];

  Print["  Valid (strict interior): ",nValid," / ",nCP,
        "  (total: ",tAll," s)"];
  Print["  -> ", If[nValid==nCP,"*** ALL VALID ***",
                   If[nValid==0,"*** ALL INVALID ***",
                     "Mixed: "<>ToString[nValid]<>"/"<>
                     ToString[nCP]]]];

  If[nCP<=24,
    Print["  Per-root (V=valid,X=invalid): ",
      Map[If[#[[4]],"V","X"]&, results]]];

  {name, nCP, nValid, dualNote}
];

(* ================================================================
   Main: all 6 regular 4-polytopes
   Note: 120-cell and 600-cell may take several minutes.
         For 600-cell, RegionMeasure on 600 hulls is expensive;
         bbox=0 roots are certified valid immediately (no hull build).
   ================================================================ *)
Print["Face-rotation BFS net check v2 (strict interior intersection)"];
Print["All 6 regular 4-polytopes\n"];

summary = {
  runPolytopeV2["5-cell",   "f5.m",   "self-dual"],
  runPolytopeV2["8-cell",   "f8.m",   "dual of 16-cell"],
  runPolytopeV2["16-cell",  "f16.m",  "dual of 8-cell"],
  runPolytopeV2["24-cell",  "f24.m",  "self-dual"],
  runPolytopeV2["120-cell", "f120.m", "dual of 600-cell"],
  runPolytopeV2["600-cell", "f600.m", "dual of 120-cell"]
};

(* ================================================================
   Final summary
   ================================================================ *)
Print["\n", StringRepeat["=",65]];
Print["SUMMARY v2 (strict interior — boundary contact ignored)"];
Print[StringRepeat["=",65]];
Print[" polytope  | cells | valid | dual            | result"];
Print[StringRepeat["-",65]];
Do[
  With[{s=r},
    Print[" ",PaddedForm[s[[1]],10],"|",PaddedForm[s[[2]],6]," |",
          PaddedForm[s[[3]],6]," | ",PaddedForm[s[[4]],16]," | ",
          Which[s[[3]]==s[[2]], "ALL VALID",
                s[[3]]==0,      "ALL INVALID",
                True, ToString[s[[3]]]<>"/"<>ToString[s[[2]]]<>" mixed"]]],
  {r,summary}];
Print[StringRepeat["=",65]];

Print["\nDual pair comparison:"];
Print["  8-cell  valid=",summary[[2,3]]," vs 16-cell valid=",summary[[3,3]]];
Print["  120-cell valid=",summary[[5,3]]," vs 600-cell valid=",summary[[6,3]]];
