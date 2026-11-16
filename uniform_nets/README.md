# Face-rotation BFS nets of the convex uniform 4-polytopes

3D nets of all **64 convex uniform 4-polytopes** (the two infinite prismatic
families excluded), produced by the face-rotation BFS unfolding described in
the accompanying paper. A breadth-first spanning tree of the cell-adjacency
graph is built from a chosen root cell, and each child cell is rotated in 4D
about the polygonal face it shares with its parent until it lands in the
parent's 3D hyperplane. The result is a set of polyhedra in 3-space.

Nets of the six *regular* 4-polytopes are in [`../face_rotation_nets/`](../face_rotation_nets/).
This directory covers the uniform case.

> **Corrected 2026-11-16.** The verdicts first published here came from an
> intersection test with two bugs (pairs of cells adjacent in the polytope
> were exempt, and some edge directions were dropped from the separating-axis
> test). They hid most of the overlaps: 21,584 valid roots became 9,872.
> `manifest.csv`, `allroots.csv` and this file now give the re-run; the nets
> themselves (the STL and OFF geometry) were right and are unchanged. See
> [`../ERRATA.md`](../ERRATA.md).

The polytopes these nets were unfolded from — vertex coordinates and the full
face lattice, one file per polytope — are in
[`../_4DData_uniform/`](../_4DData_uniform/). Start there if you want to run
the unfolding yourself rather than look at its output.

## Contents

```
stl/            67 binary STL files
off/            67 OFF files
manifest.csv    one row per file
```

**67 files for 64 polytopes.** Each polytope appears once, unfolded from root
cell 1. Three polytopes appear a second time, from roots that were chosen
with the old, buggy verdicts (see *The exceptional cases* below).

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

## ⚠️ Sixteen of these files are self-intersecting and cannot be printed

| file | polytope | pairs overlapping |
|---|---|---:|
| `x3o4o3x_root1` | $t_{0,3}\{3,4,3\}$ runcinated 24-cell | 28 |
| `x3o4o3x_root47` | $t_{0,3}\{3,4,3\}$ runcinated 24-cell | 20 |
| `x3o4x3x_root1` | $t_{0,2,3}\{3,4,3\}$ | 24 |
| `s3s4o3o_root1` | snub 24-cell | 3 |
| `o5o3x3o_root1` | $t_{2}\{5,3,3\}$ | 124 |
| `o5o3x3x_root1` | $t_{2,3}\{5,3,3\}$ | 66 |
| `o5x3x3o_root1` | $t_{1,2}\{5,3,3\}$ | 12 |
| `o5x3o3x_root1` | $t_{1,3}\{5,3,3\}$ | 16 |
| `o5x3x3x_root1` | $t_{1,2,3}\{5,3,3\}$ | 60 |
| `x5o3x3o_root1` | $t_{0,2}\{5,3,3\}$ | 4 |
| `x5x3x3o_root1` | $t_{0,1,2}\{5,3,3\}$ | 4 |
| `x5o3x3x_root1` | $t_{0,2,3}\{5,3,3\}$ runcitruncated 600-cell | 24 |
| `x5o3x3x_root14` | $t_{0,2,3}\{5,3,3\}$ runcitruncated 600-cell | 26 |
| `x5x3o3x_root1` | $t_{0,1,3}\{5,3,3\}$ runcitruncated 120-cell | 52 |
| `x5x3o3x_root748` | $t_{0,1,3}\{5,3,3\}$ runcitruncated 120-cell | 82 |
| `x5x3x3x_root1` | $t_{0,1,2,3}\{5,3,3\}$ omnitruncated 120-cell | 12 |

Their cells genuinely interpenetrate, so the meshes are **non-manifold**.
A slicer will either refuse them or produce something arbitrary in the
overlapping region. The other 51 files are valid nets. In every overlap we
inspected, the two cells are adjacent in the polytope: they lie on the two
sides of a ridge that the tree cuts.

## The exceptional cases

With the corrected test, 51 of the 64 polytopes give a valid net from every
root cell, 5 from some roots only, and 8 from none:

```
none     x3o4o3x  x3o4x3x  o5o3x3o  o5o3x3x  o5x3x3x  x5o3x3x  x5x3o3x  x5x3x3x
some     s3s4o3o    48 / 144     o5x3x3o   120 / 720    o5x3o3x     1 / 1440
         x5o3x3o  1320 / 1920    x5x3x3o  1320 / 1920
```

The extra files `x3o4o3x_root47`, `x5o3x3x_root14` and `x5x3o3x_root748` were
chosen with the old verdicts (the best root of the runcinated 24-cell, and a
failing root of each of the other two). Under the corrected test they are
simply three more failures.

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
degree. The corrected separating-axis test reproduces, tree by tree, the
overlap counts of an exact certifier (rational and $\mathbb{Q}(\sqrt5)$
arithmetic, `src/certify/`) on 4,448 BFS trees of the 120-cell and 1,200 of
the 600-cell. That certifier handles the regular polytopes only; the uniform
verdicts are floating point, with overlaps 0.006 to 0.45 deep and of positive
volume where we sampled them.

## allroots.csv

`manifest.csv` describes only the 67 nets exported here. `allroots.csv` has the
full result set: **one row per root cell of every polytope, 24,487 rows**, of
which 9,872 are valid and 14,615 self-intersect (corrected 2026-11-16; the
first version said 21,584 and 2,903).

| column | meaning |
|---|---|
| `coxeter_symbol`, `name`, `family`, `cells` | as in `manifest.csv` |
| `root` | the root cell, from 1 to `cells` |
| `verdict` | `valid` or `self-intersecting` |
| `intersecting_pairs` | pairs of cells whose interiors overlap |
| `intersection_test` | `volume` or `sat`, see below |

Thirteen polytopes have a row that is not `valid`; they are listed under
*The exceptional cases* above.

`intersection_test` is `sat` for every row: all rows were recomputed with the
corrected separating-axis test (`src/uniform/recheck_uniform_allroots.wls`,
results in `data/uniform_recheck/`).

The per-root results, with the number of candidate pairs that the CSV omits,
are in `data/uniform_recheck/<polytope>.json`.
