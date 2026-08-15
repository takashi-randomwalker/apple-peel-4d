# Combinatorial data of the 64 convex uniform 4-polytopes

One Mathematica file per polytope, holding the vertex coordinates and the full
face lattice. This is the input to the face-rotation BFS unfolding; the
resulting nets are in [`../uniform_nets/`](../uniform_nets/).

`../_4DData/` holds the same kind of data for the six *regular* 4-polytopes.
This directory covers the uniform case, in a format that differs from it in
two deliberate ways — see *Differences from `_4DData/`* below.

## Format

Each file is a list of five elements, readable with `Get`:

```mathematica
{verts, edges, faces, cellsByFace, cellAdj} = Get["x5x3x3x.m"];
```

| | |
|---|---|
| `[[1]] verts` | vertex coordinates in $\mathbb{R}^4$, machine precision |
| `[[2]] edges` | vertex-index pairs, each edge once |
| `[[3]] faces` | 2-faces as vertex-index lists, in **cyclic order** |
| `[[4]] cellsByFace` | cells as face-index lists |
| `[[5]] cellAdj` | for each cell, the cells sharing a 2-face with it |

All indices are 1-based, as Mathematica expects.

Faces are stored in cyclic order rather than sorted, because the unfolding
builds the rotation plane from the first two vertices of a shared face. Sorted
order would put antipodal vertices first for a square face and the plane would
degenerate.

Every polytope is generated with **edge length 2** and is centered at the
origin, so different polytopes are directly comparable in scale.

To get the vertices of cell `i`:

```mathematica
cellVerts[i_] := verts[[Union @@ faces[[cellsByFace[[i]]]]]]
```

## Differences from `_4DData/`

`_4DData/f*.m` has six elements:
`{verts, edges, vertAdj, faces, faceAdj, cellsByFace}`.
Here `vertAdj` and `faceAdj` are dropped and `cellAdj` is added.

**Why the two were dropped.** `face_rotation_net_all4D_v2.m` reads exactly
three fields — `raw[[1]]`, `raw[[4]]`, `raw[[6]]` — so vertex adjacency and
face adjacency are never used by anything in this repository. Across the 45
Wythoffian members they would cost 2.2 MB and 6.4 MB respectively. Both are one
line to recover if you want them:

```mathematica
vertAdj = Lookup[GroupBy[Join[edges, Reverse /@ edges], First -> Last],
                 Range[Length[verts]], {}];
faceAdj = With[{em = GroupBy[Flatten[MapIndexed[Function[{f, k},
             {Sort[#], First[k]} & /@ Partition[f, 2, 1, 1]], faces], 1],
             First -> Last]},
  Table[DeleteCases[Union @@ Lookup[em, Sort /@ Partition[faces[[k]], 2, 1, 1]],
                    k], {k, Length[faces]}]];
```

**Why cell adjacency was added.** It is the one adjacency the unfolding needs,
and `_4DData/` does not carry it: the pipeline recomputes it on every load from
`cellsByFace`. It costs 1.3 MB across the 45 Wythoffian members, a fifth of
what face adjacency costs, and it is the field a reader working on nets will
reach for first. It is equivalent to

```mathematica
cellAdj === Table[Select[Range[Length[cellsByFace]],
   # =!= i && Intersection[cellsByFace[[i]], cellsByFace[[#]]] =!= {} &],
   {i, Length[cellsByFace]}]
```

which is checked on load in `export_uniform_data.m`.

**Coordinates are machine precision, not the 30 digits the generator works in.**
The unfolding applies `N[...]` immediately on load, so the extra digits are
never used downstream and would only suggest a precision that is not there.

## index.csv

| column | meaning |
|---|---|
| `file` | file name in this directory |
| `coxeter_symbol` | ringed Coxeter diagram, e.g. `x5x3o3x`; prisms and the two diminished 600-cells are named in words |
| `name` | Schläfli t-notation, e.g. `t_{0,1,3}{5,3,3}`, or the common name where there is an unambiguous one |
| `family` | `Wythoffian`, `prism`, or `diminished 600-cell` |
| `vertices`, `edges`, `faces`, `cells` | counts |
| `mean_cell_degree` | mean over cells of the number of neighbors |
| `size_MB` | file size |

## How these were generated

| script | what it makes |
|---|---|
| `../uniform4D_wythoff.m` | the 45 Wythoffian members, from Coxeter diagrams |
| `../uniform4D_special.m` | snub 24-cell and grand antiprism, as diminishings of the 600-cell |
| `../uniform4D_prisms.m` | the 17 prisms over the Platonic and Archimedean solids |
| `../export_uniform_data.m` | this directory |

Every polytope is checked against the group-theoretic predictions for the
vertex and cell counts, the Euler relation $V-E+F-C=0$, the requirement that
each ridge lie in exactly two cells, and vertex-transitivity of the vertex
degree. The files here are re-read after writing and re-checked.

Each of the 45 Wythoffian members is determined by seven small integers — the
three Coxeter marks and the four ring bits — so this directory is a
convenience, not the source of truth. Regenerating all of it takes about
8 minutes.
