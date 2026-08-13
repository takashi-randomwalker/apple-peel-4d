# Face-rotation BFS nets of the convex uniform 4-polytopes

3D nets of all **64 convex uniform 4-polytopes** (the two infinite prismatic
families excluded), produced by the face-rotation BFS unfolding described in
the accompanying paper. A breadth-first spanning tree of the cell-adjacency
graph is built from a chosen root cell, and each child cell is rotated in 4D
about the polygonal face it shares with its parent until it lands in the
parent's 3D hyperplane. The result is a set of polyhedra in 3-space.

Nets of the six *regular* 4-polytopes are in [`../face_rotation_nets/`](../face_rotation_nets/).
This directory covers the uniform case.

## Contents

```
stl/            67 binary STL files
off/            67 OFF files
manifest.csv    one row per file
```

**67 files for 64 polytopes.** Each polytope appears once, unfolded from root
cell 1. Three polytopes appear a second time, from a root chosen to show
behaviour that root 1 does not — see *The exceptional cases* below.

## Formats

| | preserves | use |
|---|---|---|
| `stl/` | geometry only | viewing, printing |
| `off/` | polygonal faces **and** the cell decomposition | analysis, re-derivation |

**STL** is binary (not ASCII), and each cell is written as its own closed
triangulated shell. Faces are fan-triangulated, so the polygonal face
structure — which is what the unfolding is *about* — is destroyed. Use STL to
look at a net or print it, not to study it.

**OFF** keeps every face as a polygon and deduplicates vertices by position.
OFF has no notion of a 3-dimensional cell, so the cell decomposition is
recorded after the face list as comment lines

```
# cell 7: 41 42 43 44 45 46
```

giving the 0-based face indices belonging to each cell. Standard OFF viewers
ignore these lines; a parser that reads them recovers the full net.

## manifest.csv

| column | meaning |
|---|---|
| `file` | base name, without extension; the same name is used in `stl/` and `off/` |
| `coxeter_symbol` | ringed Coxeter diagram, e.g. `x5x3o3x` |
| `name` | Schläfli t-notation, e.g. `t_{0,1,3}{5,3,3}`; prisms and the two diminished 600-cells are named in words |
| `family` | `Wythoffian`, `prism`, or `diminished 600-cell` |
| `cells` | number of cells of the polytope, hence the number of possible roots |
| `root` | which cell was used as the BFS root |
| `verdict` | `valid` or `self-intersecting` |
| `intersecting_pairs` | number of pairs of cells whose interiors overlap; 0 iff valid |
| `net_vertices`, `net_faces` | size of the net as written to OFF |
| `stl_triangles` | after triangulation |

The names are given as t-notation rather than in English on purpose. Assigning
English names to all 45 Wythoffian members invites transcription errors; the
t-notation is read off the ring pattern mechanically and is unambiguous.

## ⚠️ Four of these files are self-intersecting and cannot be printed

| file | polytope | pairs overlapping |
|---|---|---|
| `x3o4o3x_root1` | runcinated 24-cell | 5 |
| `x3o4o3x_root47` | runcinated 24-cell | 3 |
| `x5o3x3x_root14` | runcitruncated 600-cell | 1 |
| `x5x3o3x_root748` | runcitruncated 120-cell | 1 |

Their cells genuinely interpenetrate, so the meshes are **non-manifold**.
A slicer will either refuse them or produce something arbitrary in the
overlapping region. They are included precisely because they are the
counterexamples: they document where the unfolding fails. The other 63 files
are valid nets with no interpenetration.

## The exceptional cases

Of the 64 polytopes, 61 give a valid net from *every* root cell. The other
three are why the extra files exist.

**Runcinated 24-cell**, `x3o4o3x` = $t_{0,3}\{3,4,3\}$, 240 cells.
No root works. `root1` is a typical failure (5 overlapping pairs); `root47`
is the best root it has, and even there 3 pairs overlap.

**Runcitruncated 600-cell**, `x5o3x3x` = $t_{0,2,3}\{5,3,3\}$, 2640 cells.
Only 47 of its 2640 roots give a valid net. `root1` happens to be one of the
47; `root14` is a failure. Note what this means: a single root tells you
almost nothing about this polytope.

**Runcitruncated 120-cell**, `x5x3o3x` = $t_{0,1,3}\{5,3,3\}$, 2640 cells.
2570 of 2640 roots are valid. `root1` is valid, `root748` is one of the 70
failures, each of which has exactly one overlapping pair.

## Scale and orientation

Every polytope is generated with edge length 2, and the net is placed with the
root cell's centroid on the $+w$ axis before projection, so nets of different
polytopes are directly comparable in size. No further scaling is applied;
scale to taste before printing.

## Reproducing

The generators and the unfolding code are in the repository root:

| script | what it makes |
|---|---|
| `uniform4D_wythoff.m` | the 45 Wythoffian members, from Coxeter diagrams |
| `uniform4D_special.m` | snub 24-cell and grand antiprism, as diminishings of the 600-cell |
| `uniform4D_prisms.m` | the 17 prisms over the Platonic and Archimedean solids |
| `bfs_sat.m` | the unfolding and the separating-axis intersection test |
| `export_uniform_nets.m` | this directory |

Each generator checks itself against the group-theoretic predictions for the
vertex and cell counts, the Euler relation $V-E+F-C=0$, the requirement that
every ridge lie in exactly two cells, and vertex-transitivity of the vertex
degree. The separating-axis test was checked to agree with a volume-based
intersection test — same verdicts and the same number of intersecting pairs —
on all 1,893 roots of the 30 polytopes of the $A_4$, $B_4$ and $F_4$ families.

## allroots.csv

`manifest.csv` describes only the 67 nets exported here. `allroots.csv` has the
full result set behind Section 4 of the paper: **one row per root cell of every
polytope, 24,487 rows**, of which 21,584 are valid and 2,903 self-intersect.

| column | meaning |
|---|---|
| `coxeter_symbol`, `name`, `family`, `cells` | as in `manifest.csv` |
| `root` | the root cell, from 1 to `cells` |
| `verdict` | `valid` or `self-intersecting` |
| `intersecting_pairs` | pairs of cells whose interiors overlap |
| `intersection_test` | `volume` or `sat`, see below |

Only three polytopes have a row that is not `valid`:

```
x3o4o3x      0 / 240  valid    runcinated 24-cell
x5o3x3x     47 / 2640 valid    runcitruncated 600-cell
x5x3o3x   2570 / 2640 valid    runcitruncated 120-cell
```

`intersection_test` records which of the two tests produced the row. The
volume test measures the volume of the intersection of two cells; the
separating-axis test is the fast equivalent used for the larger families. They
were checked to agree on the 1,893 rows where both were run, on the verdict and
on `intersecting_pairs`, and again on a sample of the root-dependent cases.
The column is kept so that the provenance of every row is explicit.

The underlying Mathematica results, including the bounding-box candidate counts
that the CSV omits, are in `wythoff_bfs_allroots.mx`, `h4_bfs_allroots.mx`,
`special_bfs_allroots.mx` and `prisms_bfs_allroots.mx`. Those are not in the
repository — regenerate them with the run scripts, or use `allroots.csv`.
