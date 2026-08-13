(* =====================================================================
   run_wythoff_h4.m
   Generate all four families including H4 (|W| = 14400) and save the
   result, so the expensive H4 pass is paid once.
   Output: wythoff_gen_all.mx  (symbol w4All)
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];

w4IncludeH4 = True;
{tGen, dummy} = AbsoluteTiming[Get[baseDir <> "uniform4D_wythoff.m"]];

Print["\ntotal generation time: ", Round[tGen, 0.1], " s"];

DumpSave[baseDir <> "wythoff_gen_all.mx", w4All];
Print["saved ", Length[w4All], " polytopes to wythoff_gen_all.mx"];

w4Distinct = DeleteDuplicatesBy[w4All, {#["counts"], #["hist"]} &];
Print["distinct by combinatorial signature: ", Length[w4Distinct]];
