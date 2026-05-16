(* peeling4Df : Apple-Peel Unfolding Algorithm for 4D Polytopes.
   4D extension of peeling3Df.
   Arguments:
     origVers : vertex coordinates {{x,y,z,w}, ...}
     faces    : faces as vertex index lists {{v1,v2,...}, ...}
     cells    : cells as FACE index lists {{f1,f2,...}, ...}
     top      : index of starting cell F1
   Returns: for each F2 candidate, {rotatedVers, order, successBool}. *)

peeling4Df[origVers_, faces_, cells_, top_] :=
  Module[
    {wHat, cc1, localVers, cCenters,
     neighbors, neighborsInit,
     order, cansInit, chosenFrom1,
     lastChosen, newChoice, cans,
     angle4D, rCans, cellVerts},

    wHat = {0, 0, 0, 1};

    (* cellVerts: vertex indices of cell i via face lookup *)
    cellVerts[i_] := Union @@ faces[[ cells[[i]] ]];

    (* Step 1: translate centroid to origin *)
    localVers = (# - Mean[origVers]) & /@ origVers;

    (* Step 2: centroid of starting cell C1 *)
    cc1 = Mean[ localVers[[ cellVerts[top] ]] ];

    (* Step 3: rotate so cc1 aligns with +w axis *)
    localVers = Which[
      Abs[(cc1 / Norm[cc1]) . wHat - 1] < 0.0001,
        localVers,
      Abs[(cc1 / Norm[cc1]) . wHat + 1] < 0.0001,
        -localVers,
      True,
        angle4D = ArcCos[Clip[cc1 . wHat / Norm[cc1], {-1, 1}]];
        RotationMatrix[angle4D, {cc1 / Norm[cc1], wHat}] . # & /@ localVers
    ];

    (* Step 4: cell centroids in rotated frame *)
    cCenters = Mean[ localVers[[ cellVerts[#] ]] ] & /@ Range[Length[cells]];

    (* Step 5: adjacency -- cells sharing >= 1 face index *)
    neighbors = Table[
      Select[
        Range[Length[cells]],
        Function[j,
          j != i && Length[Intersection[cells[[i]], cells[[j]]]] >= 1]
      ],
      {i, Length[cells]}
    ];

    (* Step 6: remove F1 from neighbor lists *)
    neighbors = neighborsInit = Complement[#, {top}] & /@ neighbors;

    (* Step 7: F2 candidates = neighbors of F1 *)
    cansInit = neighbors[[top]];

    (* Step 8: peeling loop for each F2 *)
    Table[
      order     = {top};
      neighbors = neighborsInit;
      chosenFrom1 = cansInit[[i]];
      AppendTo[order, chosenFrom1];
      neighbors = Complement[#, {chosenFrom1}] & /@ neighbors;
      lastChosen = chosenFrom1;
      newChoice  = lastChosen;

      While[
        Length[Flatten[neighbors]] != 0 &&
        Length[neighbors[[newChoice]]] != 0,

        (* left half-space: c_y*c'_x - c_x*c'_y >= 0, right-handed *)
        cans = Select[
          neighbors[[lastChosen]],
          Function[j,
            With[{c = cCenters[[lastChosen]]},
              c[[2]] * cCenters[[j]][[1]] -
              c[[1]] * cCenters[[j]][[2]] >= 0
            ]
          ]
        ];

        newChoice = Which[
          Length[cans] > 1,
            cans[[
              First[First[
                Position[cCenters[[cans, 4]], Max[cCenters[[cans, 4]]]]
              ]]
            ]],
          Length[cans] == 1,
            cans[[1]],
          True,
            rCans = neighbors[[lastChosen]];
            rCans[[
              First[First[
                Position[cCenters[[rCans, 4]], Min[cCenters[[rCans, 4]]]]
              ]]
            ]]
        ];

        lastChosen = newChoice;
        AppendTo[order, newChoice];
        neighbors = Complement[#, {newChoice}] & /@ neighbors;
      ];

      {localVers, order, Length[order] == Length[cells]}
      , {i, Length[cansInit]}
    ]

  ]
