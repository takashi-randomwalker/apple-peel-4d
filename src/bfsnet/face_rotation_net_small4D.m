(* face_rotation_net_small4D.m
   Face-rotation BFS net check for 5-cell, 8-cell, 16-cell, 24-cell.
   Self-contained. Uses polytope-generic globals (versP, facesP, cellsP).

   Reference results (already computed):
     120-cell : 120/120 VALID
     600-cell :   0/600 INVALID

   Dual pairs:
     5-cell  <-> self-dual
     8-cell  <-> 16-cell
     24-cell <-> self-dual
     120-cell<-> 600-cell  (VALID vs INVALID)
*)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";

(* ================================================================
   Polytope-independent 4D affine utilities
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
   Generic data-dependent functions
   All use globals: versP, facesP, cellsP, nCP, adjListP, adjSetP
   ================================================================ *)

(* Cell geometry *)
cVertsIdxP[i_] := Union @@ facesP[[ cellsP[[i]] ]];
cVertsP[i_]    := versP[[ cVertsIdxP[i] ]];
cCentroidP[i_] := Mean @ N @ cVertsP[i];

(* Shared face vertices (the polygon two adjacent cells share) *)
sharedFaceVertsP[i_,j_] :=
  versP[[ facesP[[ First[Intersection[cellsP[[i]], cellsP[[j]]]] ]] ]];

(* Adjacency test *)
adjQP[i_,j_] := KeyExistsQ[adjSetP, Sort[{i,j}]];

(* Face-up rotation: align root cell r's centroid to +w axis *)
faceUpForRootP[r_, locV_] := Module[
  {ccIdx, cc, wHat={0.,0.,0.,1.}, angle, R},
  ccIdx = Union @@ facesP[[ cellsP[[r]] ]];
  cc    = Mean[locV[[ccIdx]]];
  Which[
    Abs[(cc/Norm[cc]).wHat - 1] < 0.0001, locV,
    Abs[(cc/Norm[cc]).wHat + 1] < 0.0001, -locV,
    True,
      angle = ArcCos[Clip[cc.wHat/Norm[cc], {-1.,1.}]];
      R     = RotationMatrix[angle, {cc/Norm[cc], wHat}];
      R.#&/@locV
  ]
];

(* Per-edge unfolding transform (works for any shared face shape) *)
unfoldTFP[pIdx_,cIdx_,pTF_] := Module[
  {fv,fvT,cf,u1,u2,rel,nPar,nChd,v,perp,e1,e2,cosT,theta},
  fv   = N @ sharedFaceVertsP[pIdx, cIdx];
  fvT  = applyAff[pTF,#]&/@fv;
  cf   = Mean[fvT];
  rel  = (#-cf)&/@fvT;
  u1   = Normalize[rel[[1]]];
  u2   = Normalize[rel[[2]] - (rel[[2]].u1)u1];
  v    = applyAff[pTF, N@cCentroidP[pIdx]] - cf;
  nPar = Normalize[v - (v.u1)u1 - (v.u2)u2];
  v    = applyAff[pTF, N@cCentroidP[cIdx]] - cf;
  nChd = Normalize[v - (v.u1)u1 - (v.u2)u2];
  cosT  = nChd.(-nPar);
  e1    = nChd;
  perp  = -nPar - cosT*nChd;
  e2    = If[Norm[perp]>10^-8, Normalize[perp],
             Normalize@First@NullSpace[{u1,u2,nChd}]];
  theta = ArcCos[Clip[cosT,{-1,1}]];
  composeAff[affRotAroundPt[cf,e1,e2,theta], pTF]
];

(* BFS unfolding from root *)
bfsUnfoldP[root_] := Module[
  {q={root}, vis={root}, tfs=Association[root->idAff4], par,nb},
  While[q=!={},
    par=First[q]; q=Rest[q];
    Do[nb=neigh;
       If[!MemberQ[vis,nb],
         AppendTo[vis,nb]; AppendTo[q,nb];
         tfs[nb]=unfoldTFP[par,nb,tfs[par]]],
       {neigh,adjListP[[par]]}]];
  tfs
];

(* Full check for one root: face-up -> BFS -> bbox -> hull intersection *)
checkRootP[r_] := Module[
  {tfsR, unf4D, unf3D, wSpread,
   bMinR, bMaxR, bboxList, nBBox, hulls, trueHits},

  versP   = faceUpForRootP[r, versCenteredP];  (* update global versP *)
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

  If[nBBox==0,
    Return[{r,0,0,True}]];

  hulls     = Table[ConvexHullMesh[unf3D[[i]]], {i,nCP}];
  trueHits  = Select[bboxList, Function[p,
    With[{i=p[[1]],j=p[[2]]},
      AnyTrue[unf3D[[j]], RegionMember[hulls[[i]],#]&] ||
      AnyTrue[unf3D[[i]], RegionMember[hulls[[j]],#]&]]]];

  {r, nBBox, Length[trueHits], Length[trueHits]==0}
];

(* ================================================================
   Run one polytope: load, build adjacency, check all roots
   ================================================================ *)
runPolytope[name_, dataFile_, dualNote_] := Module[
  {raw, tLoad, tAdj, tAll, results, nValid},

  Print["\n", StringRepeat["-",55]];
  Print["Polytope: ", name, "  (", dualNote, ")"];
  Print[StringRepeat["-",55]];

  (* Load *)
  {tLoad, raw} = AbsoluteTiming[Get[baseDir<>"_4DData/"<>dataFile]];
  versP  = raw[[1]]; facesP = raw[[4]]; cellsP = raw[[6]];
  nCP    = Length[cellsP];
  versCenteredP = (#-Mean[N@versP])&/@N@versP;
  Print["  vertices=",Length[versP],"  faces=",Length[facesP],
        "  cells=",nCP,"  (loaded in ",tLoad," s)"];

  (* Adjacency *)
  {tAdj, adjListP} = AbsoluteTiming[
    Table[Select[Range[nCP],
      Function[j, j=!=i &&
        Length[Intersection[cellsP[[i]],cellsP[[j]]]]>=1]],
      {i,nCP}]];
  adjSetP = Association@Flatten[
    Table[Sort[{i,j}]->True,{i,nCP},{j,adjListP[[i]]}],1];
  Print["  adj degree: min=",Min[Length/@adjListP],
        " max=",Max[Length/@adjListP],
        "  (built in ",tAdj," s)"];

  (* Check all roots *)
  {tAll, results} = AbsoluteTiming[
    Table[checkRootP[r], {r,1,nCP}]];
  nValid = Length[Select[results, #[[4]]&]];

  Print["  Valid: ",nValid," / ",nCP,
        "  (", tAll," s total)"];
  Print["  -> ",
    If[nValid==nCP, "*** ALL VALID ***",
      If[nValid==0,  "*** ALL INVALID ***",
        "Mixed: "<>ToString[nValid]<>"/"<>ToString[nCP]]]];

  (* Show first few per-root results for small polytopes *)
  If[nCP <= 24,
    Print["  Per-root: ",
      Map[{#[[1]], If[#[[4]],"V","X"]}&, results]]];

  {name, nCP, nValid, dualNote}
];

(* ================================================================
   Main: check all four polytopes
   ================================================================ *)
Print["Face-rotation BFS net check: 5 / 8 / 16 / 24-cell"];
Print["(120-cell=120/120 VALID,  600-cell=0/600 INVALID for reference)\n"];

summary = {
  runPolytope["5-cell",  "f5.m",  "self-dual"],
  runPolytope["8-cell",  "f8.m",  "dual of 16-cell"],
  runPolytope["16-cell", "f16.m", "dual of 8-cell"],
  runPolytope["24-cell", "f24.m", "self-dual"]
};

(* ================================================================
   Summary table
   ================================================================ *)
Print["\n", StringRepeat["=",60]];
Print["SUMMARY: Face-rotation BFS net validity (all regular 4-polytopes)"];
Print[StringRepeat["=",60]];
Print[" polytope  | cells | valid | dual pair        | result"];
Print[StringRepeat["-",60]];
Do[
  With[{s=r},
    Print[" ",PaddedForm[s[[1]],10],"|",PaddedForm[s[[2]],6]," |",
          PaddedForm[s[[3]],6]," | ",PaddedForm[s[[4]],17]," | ",
          If[s[[3]]==s[[2]],"ALL VALID",
            If[s[[3]]==0,"ALL INVALID",
              ToString[s[[3]]]<>"/"<>ToString[s[[2]]]<>" mixed"]]]],
  {r,summary}];
Print[" 120-cell  |   120 |   120 | dual of 600-cell | ALL VALID"];
Print[" 600-cell  |   600 |     0 | dual of 120-cell | ALL INVALID"];
Print[StringRepeat["=",60]];
Print["\nDuality hypothesis check:"];
Print["  8-cell  valid=",summary[[2,3]],
      "  16-cell valid=",summary[[3,3]],
      If[summary[[2,3]]!=summary[[3,3]],
        "  -> results DIFFER (dual pair shows contrast)",
        "  -> results same"]];
Print["  120-cell: 120/120  600-cell: 0/600  -> dual pair CONTRASTS"];
