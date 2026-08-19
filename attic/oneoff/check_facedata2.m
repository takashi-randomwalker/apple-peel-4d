(* adj の実際の構造を直接確認 *)
scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "peeling3DLoxo.m"];
Get[scriptDir <> "archimedean_faceup_results.mx"];

name = "TruncatedTetrahedron";
{vers, faces} = p3lExtractPolyhedron[name];
faceType = Length /@ faces;

Print["face types (faceType): ", faceType];
Print["Total faces: ", Length[faces]];

(* allResults の details から F1,F2 の実際の型を確認 *)
details = allResults[name]["details"]["R1"];
Print["\nFirst 10 pairs (top, f2, success) with face types:"];
Do[
  {top, f2, ord, suc} = details[[i]];
  Print["  i=", i,
        "  top=", top, "(type=", faceType[[top]], ")",
        "  f2=", f2, "(type=", faceType[[f2]], ")",
        "  success=", suc]
, {i, Min[10, Length[details]]}];

(* ペアタイプの全集計 *)
pairTypes = {faceType[[#[[1]]]], faceType[[#[[2]]]]} & /@ details;
Print["\nAll (type(F1),type(F2)) pair counts: ", Tally[Sort[pairTypes]]];
Print["Total pairs: ", Length[details]];
