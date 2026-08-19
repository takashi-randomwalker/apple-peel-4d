(* debug_120cell_timing.m — time each step for the 120-cell overlap check *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];

Get[scriptDir <> "peeling4Df4.m"];
Get[scriptDir <> "unfold3DExport.m"];

dataDir = scriptDir <> "_4DData/";

raw   = Get[dataDir <> "f120.m"];
vers  = N[raw[[1]]];
faces = raw[[4]];
cells = raw[[6]];
nc    = Length[cells];
Print["120-cell: ", nc, " cells, ", Length[vers], " vertices, ", Length[faces], " faces"];

(* time a single unfoldCellLocal3D call *)
t0 = AbsoluteTime[];
emb1 = unfoldCellLocal3D[1, vers, cells, faces];
Print["unfoldCellLocal3D[1]: ", Round[AbsoluteTime[] - t0, 0.001], " s"];

t0 = AbsoluteTime[];
emb2 = unfoldCellLocal3D[2, vers, cells, faces];
Print["unfoldCellLocal3D[2]: ", Round[AbsoluteTime[] - t0, 0.001], " s"];

(* pre-compute all 120 local3D *)
t0 = AbsoluteTime[];
cellLocal3Ds = Table[unfoldCellLocal3D[ci, vers, cells, faces], {ci, nc}];
Print["All ", nc, " cellLocal3D: ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* run one peeling to get an order *)
t0 = AbsoluteTime[];
ans1 = peeling4Df4[vers, faces, cells, 1, "RZ"];
Print["peeling4Df4 C1=1: ", Round[AbsoluteTime[] - t0, 0.01], " s"];

success = Select[ans1, Last[#] &];
Print["Successful C2 count: ", Length[success]];
If[Length[success] == 0, Print["No success, exiting"]; Exit[]];

order = success[[1, 2]];
Print["First order length: ", Length[order]];

(* time unfoldTo3D *)
t0 = AbsoluteTime[];
cellEmbs = unfoldTo3D[vers, faces, cells, order];
Print["unfoldTo3D (original): ", Round[AbsoluteTime[] - t0, 0.01], " s"];

t0 = AbsoluteTime[];
cellEmbsFast = unfoldTo3DFast[cells, faces, order, cellLocal3Ds];
Print["unfoldTo3DFast (pre-computed): ", Round[AbsoluteTime[] - t0, 0.01], " s"];

(* time unfoldCellSATData for one cell *)
t0 = AbsoluteTime[];
satd1 = unfoldCellSATData[order[[1]], cells, faces, cellEmbs[[1]]];
Print["unfoldCellSATData (1 cell): ", Round[AbsoluteTime[] - t0, 0.001], " s"];

(* time building all SAT data *)
t0 = AbsoluteTime[];
cds = Table[unfoldCellSATData[order[[k]], cells, faces, cellEmbs[[k]]], {k, nc}];
Print["All ", nc, " unfoldCellSATData: ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* time checkOverlaps *)
t0 = AbsoluteTime[];
ovs = checkOverlaps[order, cells, faces, cellEmbs];
Print["checkOverlaps: ", Round[AbsoluteTime[] - t0, 0.1], " s  (overlaps: ", Length[ovs], ")"];

Print["Done."];
