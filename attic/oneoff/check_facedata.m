(* face representation の診断 *)
scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "peeling3DLoxo.m"];

diag[name_] := Module[{vers, faces, faceType, adj, typeCounts, adjTypes},
  {vers, faces} = p3lExtractPolyhedron[name];
  faceType = Length /@ faces;
  typeCounts = Counts[faceType];
  adj = p3lBuildAdj[faces];
  (* adj[[f]] = face f に隣接する面インデックスのリスト *)
  adjTypes = Table[
    Sort[{faceType[[f]], faceType[[#]]}] & /@ adj[[f]],
    {f, Length[faces]}];
  adjTypeFlat = Tally[Sort[Flatten[adjTypes, 1]]];
  Print["\n", name];
  Print["  faces: ", Length[faces], "  face-type counts: ", Normal[typeCounts]];
  Print["  undirected adj-type pairs: ", adjTypeFlat];
  (* directed pairs grouped by face type *)
  dirPairs = Flatten[Table[
    {faceType[[f]], faceType[[g]]} & /@ adj[[f]],
    {f, Length[faces]}], 1];
  Print["  directed adj-type pairs: ", Tally[Sort[dirPairs]]]
];

diag["TruncatedTetrahedron"];
diag["SnubCube"];
diag["TruncatedCuboctahedron"];
