# Apple-Peel Unfolding in Three and Four Dimensions

Mathematica implementation and 3D-printable nets for the Apple-Peel Unfolding algorithm applied to Platonic solids, Archimedean solids, and regular convex 4-polytopes.

## Overview

Apple-Peel Unfolding is a greedy algorithm that selects faces (or cells) of a polyhedron (or 4-polytope) one at a time in a spiral order, producing a net analogous to peeling an apple in a single continuous strip.

Two selection rules are implemented:

- **RS (Spiral rule)**: minimum signed determinant (sharpest clockwise turn)
- **RZ (Zonal rule)**: maximum coordinate along the peeling axis

### Results Summary (RZ rule)

| Object | Pairs | Result |
|--------|------:|--------|
| 5 Platonic solids | — | Perfect (100%) |
| 13 Archimedean solids | — | 3 Perfect, 3 Possible, 7 Impossible |
| 5-cell | 20 | **Perfect** |
| 8-cell | 48 | **Perfect** |
| 16-cell | 64 | Possible (20/64) |
| 24-cell | 192 | **Perfect** |
| 120-cell | 1,440 | **Perfect** |
| 600-cell | 2,400 | Impossible |

## Requirements

- Mathematica 13 or later

## File Structure

```
peeling3DLoxo.m      — 3D algorithm (RS and RZ rules)
peeling4Df4.m        — 4D algorithm (global c1–c2 reference)
unfold3DExport.m     — 2D/3D net construction and STL export
_4DData/             — Regular 4-polytope data (f5.m … f600.m)
STLs_v4/             — 280 valid 3D-printable nets (no self-intersection)
run_*.m              — Batch computation scripts
visualize_*.m        — Visualization scripts
```

## Usage

### 3D Polyhedra

```mathematica
Get["peeling3DLoxo.m"];
{vers, faces} = p3lExtractPolyhedron["Dodecahedron"];
result = p3lRunAll[vers, faces, "maxz"];  (* RZ rule *)
```

### 4D Polytopes

```mathematica
Get["peeling4Df4.m"];
Get["_4DData/f120.m"];  (* 120-cell *)
result = p4fRunAll[vers, faces, cells, "RZ"];
```

### 3D-Printable Nets

The `STLs_v4/` directory contains 280 STL files — all valid (non-self-intersecting) nets for the 5-cell, 8-cell, 16-cell, and 24-cell under the RZ rule.

File naming: `unfold_{polytope}_C1_{i}_C2_{j}.stl`

## Reference

T. Yoshino and S. Chaidee,
*Apple-Peel Unfolding in Three and Four Dimensions: Spiral and Zonal Selection Rules*,
preprint, 2026.

Companion paper (3D Archimedean and Catalan solids):
T. Yoshino and S. Chaidee, arXiv:2604.16204, 2026.

## License

MIT License — see [LICENSE](LICENSE) for details.
