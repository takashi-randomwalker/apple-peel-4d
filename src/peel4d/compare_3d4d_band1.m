(* ================================================================
   compare_3d4d_band1.m
   120-cell の第1バンド訪問順序と正十二面体3D皮むき順序の比較
   ================================================================ *)

$dir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
dataDir = $dir <> "_4DData/";

(* --- 120-cell データ --- *)
raw = Get[dataDir <> "f120.m"];
vers = N[raw[[1]]]; faces = raw[[4]]; cells = raw[[6]];
nc = Length[cells];
Print["120-cell: nc = ", nc];

(* --- 4D peeling 読み込み --- *)
Get[$dir <> "peeling4Df4.m"];

(* --- セル重心・隣接 --- *)
cv[i_] := Union @@ faces[[ cells[[i]] ]];
adj = Table[Select[Range[nc], Function[j, j != i &&
  Length[Intersection[cells[[i]], cells[[j]]]] >= 1]], {i, nc}];

top = 1;
band8 = adj[[top]];
Print["Band-8 size: ", Length[band8], "  cells: ", Sort[band8]];

(* --- C1=1 で RZ 実行（12 C2 候補分） --- *)
Print["\nRunning peeling4Df4 for C1=1 (RZ)..."];
res = peeling4Df4[vers, faces, cells, top, "RZ"];
Print["Done. ", Length[res], " pairs."];

(* --- 検証1: 全 C2 で k=1..13 が {C1} ∪ band8 か --- *)
Print["\n=== Verification: first 13 steps = C1 + band-8? ==="];
allOK = True;
Do[
  ord = res[[i]][[2]];
  match = Sort[ord[[1;;13]]] === Sort[Join[{top}, band8]];
  If[!match, Print["C2=", ord[[2]], " MISMATCH"]; allOK = False],
  {i, Length[res]}
];
If[allOK, Print["ALL 12 C2 choices: first 13 steps = C1 + band-8. OK"]];

(* --- 検証2: band-8 の w 値（全員等値のはず） --- *)
localVers1 = res[[1]][[1]];
cC = Mean[localVers1[[ cv[#] ]]] & /@ Range[nc];
c1v = cC[[top]];
Print["\n=== Band-8 w-coordinates (should all be equal) ==="];
wvals = cC[[band8, 4]];
Print["w values: ", NumberForm[#, {6,5}] & /@ N[wvals]];
Print["spread:   ", N[Max[wvals] - Min[wvals], 4]];

(* --- 検証3: band-8 の xy角度と訪問順序（C2別） --- *)
Print["\n=== Band-8 visit order vs xy angle for each C2 ==="];
Do[
  ord = res[[i]][[2]];
  c2 = ord[[2]];
  visitSeq = ord[[2;;13]]; (* band-8 cells in visit order *)
  angles = N[ArcTan[cC[[#,1]], cC[[#,2]]] * 180/Pi, 4] & /@ visitSeq;
  Print["C2=", c2, ": visit order angles (deg) = ",
        NumberForm[#, {6,2}] & /@ angles],
  {i, Length[res]}
];

(* --- 3D 比較: 正十二面体の面重心 z 値と RZ 順序 --- *)
Print["\n=== 3D Dodecahedron face centroids z-values ==="];
(* Mathematica 組み込み十二面体 *)
dVers = N[PolyhedronData["Dodecahedron", "VertexCoordinates"]];
dFaces = PolyhedronData["Dodecahedron", "FaceIndices"];
nf = Length[dFaces];
Print["Dodecahedron: ", nf, " faces"];
fCentroids = Mean[dVers[[#]]] & /@ dFaces;

(* face-up: face 1 の重心を +z に揃える *)
f1c = fCentroids[[1]];
zHat = {0,0,1};
If[Abs[(f1c/Norm[f1c]) . zHat - 1] < 0.001,
  localDVers = dVers,
  angle3D = ArcCos[Clip[f1c . zHat / Norm[f1c], {-1,1}]];
  R3D = RotationMatrix[angle3D, {f1c/Norm[f1c], zHat}];
  localDVers = R3D . # & /@ dVers
];
fc3D = Mean[localDVers[[#]]] & /@ dFaces;

(* 3D RZ 順序: max z, 次に min Det で *)
Get[$dir <> "peeling3DLoxo.m"];
res3D = p3lRunAll[localDVers, dFaces, "RZ"];
(* 成功しているものの最初の順序を取得 *)
successIdx = Select[Range[Length[res3D]], res3D[[#]][[3]] &];
Print["3D RZ successful pairs: ", Length[successIdx], " / ", Length[res3D]];
If[Length[successIdx] > 0,
  ord3D = res3D[[successIdx[[1]]]][[2]];
  Print["3D RZ face visit order (pair ", successIdx[[1]], "): ", ord3D];
  Print["z-values in visit order: ",
    NumberForm[#, {6,4}] & /@ (fc3D[[#, 3]] & /@ ord3D)]
];

(* band-8 の w 値は全等値 → 訪問順は geo-score（Det）で決まる *)
(* 4D Det → -w1 * Det_xyz{c2,ck,cj} これは実質3D行列式 *)
Print["\n=== Structural conclusion ==="];
Print["Band-8 w-spread = ", N[Max[wvals] - Min[wvals], 2]];
Print["Since all band-8 cells share the same w, the RZ ordering"];
Print["within band-8 reduces to the geo-score (Det) criterion,"];
Print["which equals -w1 * Det_xyz{c2,ck,cj}: a 3D determinant."];

Print["\nDone."];
