"""
certify_bfs_exact.py
Exact certification that every BFS net of a regular 4-polytope is valid.

Everything is done over Q with fractions.Fraction; no floating point is
used anywhere.  Independent of the Mathematica pipeline: it reads only
the combinatorial data in _4DData/ (the 5-cell is re-embedded, see below).

Unfolding without trigonometry.  A regular cell is determined by one of
its faces together with the side on which it lies.  So when child c is
unfolded across the face f it shares with its parent p, the placed child
is the isometry of aff(c) into the root's 3-flat that fixes the placed
vertices of f and sends the centroid of c to the reflection of the placed
centroid of p in the plane of f.  Four affinely independent points fix
the map, and a reflection needs no normalized normal, so the placed
coordinates stay in Q whenever the input coordinates do.  Every placement
is checked exactly: it lies in the root's 3-flat and preserves the
squared length of every edge of the cell.

All BFS trees at once.  The placement of a cell depends only on its path
to the root.  A pair (path P to u, path Q to v) occurs together in some
BFS tree iff P and Q assign the same parent to every cell they share;
the remaining parent choices are then free.  So certifying every such
compatible pair, for every pair of cells, certifies every BFS tree --
including the pairs joined by a tree edge, which Definition 2 exempts.

Separation.  Interiors of two full-dimensional convex polytopes A, B are
disjoint iff 0 is not interior to A - B, iff some facet normal u of A - B
has max_A u.x <= min_B u.x.  Facet normals of A - B are face normals of A
or B or cross products of an edge of each.  We try the centroid
difference first (any separating u is a certificate), then that finite
set; if all fail, the interiors provably overlap.

Usage: python3 -I certify_bfs_exact.py <5|8|16|24> [root ...]
"""

import ast
import gzip
import itertools
import json
import os
import sys
import time
from fractions import Fraction as Q

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # -I drops the script dir
from qsqrt5 import PHI, R5, SQRT5

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))


# ---------------------------------------------------------------- data

def load_polytope(name):
    """Return (verts, faces, cells, edges); faces/cells/edges 0-based."""
    if name == "5":
        # The file's coordinates involve GoldenRatio.  The regular 5-cell is
        # congruent (up to scale) to conv(e_1..e_5) in R^5, which is rational.
        n = 5
        verts = [[Q(int(i == j)) for j in range(n)] for i in range(n)]
        faces = [list(t) for t in itertools.combinations(range(n), 3)]
        cells = [[k for k, f in enumerate(faces) if set(f) <= set(c)]
                 for c in itertools.combinations(range(n), 4)]
        edges = [tuple(e) for e in itertools.combinations(range(n), 2)]
        return verts, faces, cells, edges
    path = os.path.join(BASE, "_4DData", "f%s.m" % name)
    with open(path) as fh:
        txt = fh.read()
    if "GoldenRatio" in txt or "Sqrt" in txt:
        # 120- and 600-cell: coordinates are monomials in GoldenRatio and Sqrt[5]
        for m, py in [("Sqrt[5]", "SQRT5"), ("GoldenRatio^(-1)", "(PHI**-1)"),
                      ("GoldenRatio^(-2)", "(PHI**-2)"), ("GoldenRatio^2", "(PHI**2)"),
                      ("GoldenRatio", "PHI")]:
            txt = txt.replace(m, py)
        txt = txt.replace("{", "[").replace("}", "]")
        assert "^" not in txt and "Sqrt" not in txt, "unexpected expression in %s" % path
        raw = eval(txt, {"__builtins__": {}, "PHI": PHI, "SQRT5": SQRT5})
        verts = [[x if isinstance(x, R5) else R5(x) for x in v] for v in raw[0]]
    else:
        raw = ast.literal_eval(txt.replace("{", "[").replace("}", "]"))
        verts = [[Q(x) for x in v] for v in raw[0]]
    edges = sorted({tuple(sorted((a - 1, b - 1))) for a, b in raw[1]})
    faces = [[i - 1 for i in f] for f in raw[3]]
    cells = [[i - 1 for i in c] for c in raw[5]]
    return verts, faces, cells, edges


# ------------------------------------------------------- linear algebra

def sub(a, b): return [x - y for x, y in zip(a, b)]
def add(a, b): return [x + y for x, y in zip(a, b)]
def scl(s, a): return [s * x for x in a]
def dot(a, b): return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def solve(m, rhs):
    """Exact Gaussian elimination for a square nonsingular system."""
    n = len(m)
    a = [list(row) + [r] for row, r in zip(m, rhs)]
    for col in range(n):
        piv = next(i for i in range(col, n) if a[i][col] != 0)
        a[col], a[piv] = a[piv], a[col]
        for i in range(n):
            if i != col and a[i][col] != 0:
                f = a[i][col] / a[col][col]
                a[i] = [x - f * y for x, y in zip(a[i], a[col])]
    return [a[i][n] / a[i][i] for i in range(n)]


def coords_in(basis, v):
    """Coefficients of v in span(basis), and the exact residual."""
    g = [[dot(b, c) for c in basis] for b in basis]
    k = solve(g, [dot(b, v) for b in basis])
    res = v
    for kk, b in zip(k, basis):
        res = sub(res, scl(kk, b))
    return k, res


def centroid(pts):
    return scl(Q(1, len(pts)), [sum(c) for c in zip(*pts)])


def independent(pts, k):
    """Indices of k+1 affinely independent points among pts (greedy)."""
    chosen = [0]
    basis = []
    for i in range(1, len(pts)):
        d = sub(pts[i], pts[0])
        if basis:
            _, res = coords_in(basis, d)
        else:
            res = d
        if any(res):
            chosen.append(i)
            basis.append(d)
            if len(basis) == k:
                return chosen
    raise ValueError("points do not span a %d-flat" % k)


# ---------------------------------------------------------- the polytope

class Polytope:
    def __init__(self, name):
        self.name = name
        self.verts, self.faces, self.cells, edges = load_polytope(name)
        self.nC = len(self.cells)
        self.cellV = [sorted({v for f in c for v in self.faces[f]}) for c in self.cells]
        self.cellE = [[e for e in edges if e[0] in set(vs) and e[1] in set(vs)] for vs in self.cellV]
        self.cellCent = [centroid([self.verts[v] for v in vs]) for vs in self.cellV]
        faceCells = {}
        for c, fs in enumerate(self.cells):
            for f in fs:
                faceCells.setdefault(f, []).append(c)
        assert all(len(cs) == 2 for cs in faceCells.values()), "a face not in exactly 2 cells"
        self.shared = {}
        self.adj = [[] for _ in range(self.nC)]
        for f, (a, b) in faceCells.items():
            self.shared[(a, b)] = self.shared[(b, a)] = f
            self.adj[a].append(b)
            self.adj[b].append(a)
        for a in self.adj:          # ascending, as adjListP in the Mathematica code
            a.sort()
        # sanity: a regular polytope -- every edge has the same squared length
        lens = {dot(sub(self.verts[x], self.verts[y]), sub(self.verts[x], self.verts[y]))
                for x, y in edges}
        assert len(lens) == 1, "edges of different lengths"

    def queue_tree(self, r):
        """The tree T_r of the paper: queue BFS, neighbors in ascending order,
        each cell takes the first cell that discovers it as its parent
        (bfsUnfoldP in src/bfsnet/face_rotation_net_all4D_v2.m)."""
        parent = [None] * self.nC
        seen = {r}
        q = [r]
        for p in q:
            for nb in self.adj[p]:
                if nb not in seen:
                    seen.add(nb)
                    parent[nb] = p
                    q.append(nb)
        return parent

    def layers(self, r):
        dist = [-1] * self.nC
        dist[r] = 0
        q = [r]
        for p in q:
            for nb in self.adj[p]:
                if dist[nb] < 0:
                    dist[nb] = dist[p] + 1
                    q.append(nb)
        opts = [[] if i == r else sorted(nb for nb in self.adj[i] if dist[nb] == dist[i] - 1)
                for i in range(self.nC)]
        return dist, opts


# ------------------------------------------------------------ unfolding

class Unfolder:
    """Exact placements of cells along root paths, memoized by path."""

    def __init__(self, P, r):
        self.P = P
        self.r = r
        rv = [P.verts[v] for v in P.cellV[r]]
        idx = independent(rv, 3)
        self.origin = rv[idx[0]]
        self.basis = [sub(rv[i], rv[idx[0]]) for i in idx[1:]]   # directions of the root 3-flat
        self.memo = {(r,): {v: P.verts[v] for v in P.cellV[r]}}
        self.chart = {}

    def place(self, path):
        """path = (r, ..., v); returns {vertex: placed point in R^d}."""
        if path in self.memo:
            return self.memo[path]
        P = self.P
        parent = self.place(path[:-1])
        p, c = path[-2], path[-1]
        f = P.shared[(p, c)]
        fv = P.faces[f]
        # reflection of the parent's placed centroid in the plane of f, within the root 3-flat
        a = centroid(list(parent.values()))
        P0, P1, P2 = (parent[fv[i]] for i in independent([parent[v] for v in fv], 2))
        e = [sub(P1, P0), sub(P2, P0)]
        n = None
        for b in self.basis:
            _, res = coords_in(e, b)
            if any(res):
                n = res
                break
        a2 = sub(a, scl(2 * dot(sub(a, P0), n) / dot(n, n), n))
        # affine map from aff(c): three independent face vertices + centroid
        i0, i1, i2 = independent([P.verts[v] for v in fv], 2)
        src = [P.verts[fv[i0]], P.verts[fv[i1]], P.verts[fv[i2]], P.cellCent[c]]
        dst = [parent[fv[i0]], parent[fv[i1]], parent[fv[i2]], a2]
        sb = [sub(s, src[0]) for s in src[1:]]
        db = [sub(d, dst[0]) for d in dst[1:]]
        out = {}
        for v in P.cellV[c]:
            k, res = coords_in(sb, sub(P.verts[v], src[0]))
            assert not any(res), "vertex outside the cell's 3-flat"
            out[v] = add(dst[0], add(add(scl(k[0], db[0]), scl(k[1], db[1])), scl(k[2], db[2])))
        for v in fv:                                   # the hinge is fixed pointwise
            assert out[v] == parent[v], "shared face moved"
        for (x, y) in P.cellE[c]:                      # isometry on every edge
            d0 = sub(P.verts[x], P.verts[y])
            d1 = sub(out[x], out[y])
            assert dot(d0, d0) == dot(d1, d1), "edge length changed"
        for v in out.values():                         # stays in the root 3-flat
            _, res = coords_in(self.basis, sub(v, self.origin))
            assert not any(res), "left the root 3-flat"
        self.memo[path] = out
        return out

    def chart3(self, path):
        """Placed cell in rational affine coordinates of the root 3-flat."""
        if path not in self.chart:
            pl = self.place(path)
            self.chart[path] = {v: coords_in(self.basis, sub(x, self.origin))[0] for v, x in pl.items()}
        return self.chart[path]


# ----------------------------------------------------------- separation

def cell_geometry(P, c, pts):
    """Points, face normals and edge directions of a placed cell (3D chart)."""
    vs = P.cellV[c]
    X = [pts[v] for v in vs]
    normals = []
    for f in P.cells[c]:
        fp = [pts[v] for v in P.faces[f]]
        i0, i1, i2 = independent(fp, 2)
        normals.append(cross(sub(fp[i1], fp[i0]), sub(fp[i2], fp[i0])))
    edges = [sub(pts[x], pts[y]) for x, y in P.cellE[c]]
    return X, normals, edges


def separates(u, A, B):
    """+1 if u.A <= u.B, -1 if u.B <= u.A, 0 otherwise."""
    pa = [dot(u, x) for x in A]
    pb = [dot(u, x) for x in B]
    if max(pa) <= min(pb):
        return 1
    if max(pb) <= min(pa):
        return -1
    return 0


def certificate(gA, gB):
    """A rational u with max u.A <= min u.B, or None if interiors overlap."""
    A, nA, eA = gA
    B, nB, eB = gB
    cand = [sub(centroid(B), centroid(A))] + nA + nB
    cand += [w for a in eA for b in eB for w in [cross(a, b)] if any(w)]
    for u in cand:
        s = separates(u, A, B)
        if s:
            return u if s > 0 else scl(-1, u)
    return None


# ----------------------------------------------------------- the driver

def root_paths(r, dist, opts):
    """All root paths to every cell, as tuples (r, ..., v)."""
    paths = {r: [(r,)]}
    for v in sorted(range(len(dist)), key=lambda i: dist[i]):
        if v != r:
            paths[v] = [pp + (v,) for p in opts[v] for pp in paths[p]]
    return paths


def compatible(P, Q_):
    """Two root paths agree on the parent of every cell they share."""
    mp = {P[i]: P[i - 1] for i in range(1, len(P))}
    for i in range(1, len(Q_)):
        w = Q_[i]
        if w in mp and mp[w] != Q_[i - 1]:
            return False
    return True


def certify_root(Pt, r, out_cert=None):
    dist, opts = Pt.layers(r)
    U = Unfolder(Pt, r)
    paths = root_paths(r, dist, opts)
    geo = {}
    for v, ps in paths.items():
        for p in ps:
            geo[p] = cell_geometry(Pt, v, U.chart3(p))
    nTrees = 1
    for v in range(Pt.nC):
        if v != r:
            nTrees *= len(opts[v])
    if out_cert is not None:
        # a functional c on the chart is (B G^-1 c).(x - origin) on R^d
        G = [[dot(b, b2) for b2 in U.basis] for b in U.basis]

        def ambient(c):
            w = solve(G, c)
            return [sum(wi * b[k] for wi, b in zip(w, U.basis)) for k in range(len(U.origin))]

        normals = []                     # complete the root 3-flat's directions to R^d
        for k in range(len(U.origin)):
            e = [Q(int(i == k)) for i in range(len(U.origin))]
            _, res = coords_in(U.basis + normals, e)
            if any(res):
                normals.append(res)
        out_cert.update({"root": r, "flat_normals": [[str(x) for x in n] for n in normals],
                         "placements": [], "pairs": []})
        for v, ps in paths.items():
            for p in ps:
                out_cert["placements"].append(
                    [list(p), {str(w): [str(x) for x in pt] for w, pt in U.place(p).items()}])
    checked = 0
    overlaps = []
    for u, v in itertools.combinations(range(Pt.nC), 2):
        for p in paths[u]:
            for q in paths[v]:
                if not compatible(p, q):
                    continue
                cert = certificate(geo[p], geo[q])
                checked += 1
                if cert is None:
                    overlaps.append((p, q))
                elif out_cert is not None:
                    out_cert["pairs"].append([list(p), list(q), [str(x) for x in ambient(cert)]])
    nPaths = sum(len(ps) for ps in paths.values())
    return {"root": r, "trees": nTrees, "placements": nPaths,
            "pairs": checked, "overlaps": len(overlaps), "examples": overlaps[:5]}


def certify_tree(Pt, r, parent):
    """Certify one tree: every pair of cells, tree edges included.

    Sweep and prune along chart coordinate 0: once a cell's maximum is <= the
    minimum of the cell being added, it is separated -- by that coordinate,
    a valid certificate -- from it and from every later cell, and leaves the
    active list.  Remaining pairs first try the other two coordinates, then
    the full separating-axis test."""
    U = Unfolder(Pt, r)

    def path(v):
        out = [v]
        while out[-1] != r:
            out.append(parent[out[-1]])
        return tuple(reversed(out))

    charts = [U.chart3(path(v)) for v in range(Pt.nC)]
    X = [[charts[v][w] for w in Pt.cellV[v]] for v in range(Pt.nC)]
    lo = [[min(p[k] for p in X[v]) for k in range(3)] for v in range(Pt.nC)]
    hi = [[max(p[k] for p in X[v]) for k in range(3)] for v in range(Pt.nC)]
    geo = {}

    def g(v):
        if v not in geo:
            geo[v] = cell_geometry(Pt, v, charts[v])
        return geo[v]

    order = sorted(range(Pt.nC), key=lambda v: lo[v][0])
    active = []
    n_axis = n_sat = 0
    overlaps = []
    for b in order:
        active = [a for a in active if not hi[a][0] <= lo[b][0]]
        for a in active:
            if any(hi[a][k] <= lo[b][k] or hi[b][k] <= lo[a][k] for k in range(3)):
                n_axis += 1
                continue
            n_sat += 1
            if certificate(g(a), g(b)) is None:
                overlaps.append((a, b))
        active.append(b)
    n_pairs = Pt.nC * (Pt.nC - 1) // 2
    n_sweep = n_pairs - n_axis - n_sat          # dropped from the active list: coordinate 0
    assert n_sweep >= 0
    depth = max(len(path(v)) for v in range(Pt.nC)) - 1
    return {"root": r, "pairs": n_pairs, "sweep": n_sweep, "axis": n_axis, "sat": n_sat,
            "overlaps": len(overlaps), "examples": overlaps[:5], "depth": depth}


def _tree_worker(args):
    name, r, parent = args
    Pt = _tree_worker.cache.get(name) or _tree_worker.cache.setdefault(name, Polytope(name))
    t = time.time()
    if parent is None:
        parent = Pt.queue_tree(r)
    else:                                       # a given tree: check that it is a BFS tree
        _, opts = Pt.layers(r)
        assert parent[r] == -1 and all(parent[v] in opts[v] for v in range(Pt.nC) if v != r)
    res = certify_tree(Pt, r, parent)
    res["time"] = time.time() - t
    return res


_tree_worker.cache = {}


def main_tree(name, roots, jobs, trees_file=None):
    """Certify T_r, or with trees_file (gzip JSON, per root a list of parent
    arrays, -1 at the root) the BFS trees listed there."""
    from concurrent.futures import ProcessPoolExecutor
    Pt = Polytope(name)
    roots = roots or list(range(Pt.nC))
    if trees_file:
        with gzip.open(trees_file, "rt") as fh:
            given = json.load(fh)
        tasks = [(name, r, p) for r in roots for p in given[r]]
        what = "%d given BFS trees" % len(tasks)
    else:
        tasks = [(name, r, None) for r in roots]
        what = "tree T_r"
    print("%s-cell: %d cells, %s, %d roots, %d processes" % (name, Pt.nC, what, len(roots), jobs), flush=True)
    tot = {"pairs": 0, "sweep": 0, "axis": 0, "sat": 0, "overlaps": 0}
    bad = []
    t0 = time.time()
    with ProcessPoolExecutor(jobs) as ex:
        for i, res in enumerate(ex.map(_tree_worker, tasks, chunksize=4)):
            for k in tot:
                tot[k] += res[k]
            if res["overlaps"]:
                bad.append(res)
            if i < 3 or (i + 1) % max(1, len(tasks) // 20) == 0 or res["overlaps"]:
                print("  root %3d  depth %2d  pairs %7d = sweep %7d + axis %6d + SAT %5d  overlaps %d  (%.1fs)"
                      % (res["root"], res["depth"], res["pairs"], res["sweep"], res["axis"], res["sat"],
                         res["overlaps"], res["time"]), flush=True)
    print("TOTAL  roots %d  trees %d  pairs %d = sweep %d + axis %d + SAT %d  overlaps %d  invalid trees %d  (%.1fs)"
          % (len(roots), len(tasks), tot["pairs"], tot["sweep"], tot["axis"], tot["sat"], tot["overlaps"], len(bad),
             time.time() - t0))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    name = args[0]
    roots = [int(x) for x in args[1:]]
    jobs = next((int(a[7:]) for a in sys.argv if a.startswith("--jobs=")), os.cpu_count())
    trees_file = next((a[8:] for a in sys.argv if a.startswith("--trees=")), None)
    if "--tree" in sys.argv or trees_file or name in ("120", "600"):
        return main_tree(name, roots, jobs, trees_file)
    Pt = Polytope(name)
    roots = roots or list(range(Pt.nC))
    print("%s-cell: %d cells, roots %s" % (name, Pt.nC, "all" if len(args) == 1 else roots))
    tot = {"trees": 0, "pairs": 0, "overlaps": 0}
    export = next((a[9:] for a in sys.argv if a.startswith("--export=")), None)
    cert_roots = []
    t0 = time.time()
    for r in roots:
        t = time.time()
        rc = {} if export else None
        res = certify_root(Pt, r, rc)
        if export:
            rc["trees"] = res["trees"]
            cert_roots.append(rc)
        for k in tot:
            tot[k] += res[k]
        print("  root %3d  trees %8d  placements %5d  pairs %8d  overlaps %d  (%.1fs)"
              % (r, res["trees"], res["placements"], res["pairs"], res["overlaps"], time.time() - t),
              flush=True)
        if res["overlaps"]:
            print("    e.g.", res["examples"])
    print("TOTAL  trees %d  pairs certified %d  overlaps %d  (%.1fs)"
          % (tot["trees"], tot["pairs"], tot["overlaps"], time.time() - t0))
    if export:
        assert tot["overlaps"] == 0, "refusing to write a certificate with overlaps"
        doc = {"polytope": name,
               "vertices": [[str(x) for x in v] for v in Pt.verts],
               "faces": Pt.faces, "cells": Pt.cells,
               "roots": cert_roots}
        with gzip.open(export, "wt") as fh:
            json.dump(doc, fh, separators=(",", ":"))
        print("certificate written to %s (%.1f MB)" % (export, os.path.getsize(export) / 1e6))


if __name__ == "__main__":
    main()
