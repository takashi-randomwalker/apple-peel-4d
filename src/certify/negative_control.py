"""
negative_control.py
Can certify_bfs_exact.py actually see an overlap?

The 24-cell is not all-net (Seamons 2026; Devadoss-Harvey-Richter 2026):
some ridge unfolding along a chain of octahedra overlaps itself.  Such a
chain is not a BFS path, so the certifier never meets it.  Here we drop
the BFS restriction, grow simple paths of adjacent cells from a root, and
ask the same separation routine whether the last cell meets an earlier
one.  If it never reports an overlap, the certifier proves nothing.

For every overlap found we also exhibit a witness: a rational point
strictly interior to both cells (found
by averaging the exact vertices of the intersection polytope).

Usage: python3 -I negative_control.py <8|16|24> <max cells on path> [root]
"""

import itertools
import os
import sys
from fractions import Fraction as Q

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # -I drops the script dir
from certify_bfs_exact import Polytope, Unfolder, cell_geometry, certificate, centroid, cross, dot, scl, solve, sub


def strictly_inside(P, c, pts, x):
    """x is interior to the placed cell c (3D chart): every face is strict."""
    cen = centroid([pts[v] for v in P.cellV[c]])
    _, normals, _ = cell_geometry(P, c, pts)
    for f, n in zip(P.cells[c], normals):
        a = pts[P.faces[f][0]]
        s_x = dot(n, sub(x, a))
        s_c = dot(n, sub(cen, a))
        if s_x * s_c <= 0:          # on the face or on the far side
            return False
    return True


def halfspaces(P, c, pts):
    """Face inequalities n.x <= b of the placed cell, oriented outward."""
    cen = centroid([pts[v] for v in P.cellV[c]])
    _, normals, _ = cell_geometry(P, c, pts)
    hs = []
    for f, n in zip(P.cells[c], normals):
        b = dot(n, pts[P.faces[f][0]])
        if dot(n, cen) > b:
            n, b = scl(-1, n), -b
        hs.append((n, b))
    return hs


def witness(P, U, p, q):
    """Average of the vertices of A cap B: interior if A cap B is 3-dimensional."""
    A, B = U.chart3(p), U.chart3(q)
    hs = halfspaces(P, p[-1], A) + halfspaces(P, q[-1], B)
    verts = set()
    for i, j, k in itertools.combinations(range(len(hs)), 3):
        m = [hs[i][0], hs[j][0], hs[k][0]]
        if det3(m) == 0:
            continue
        x = solve(m, [hs[i][1], hs[j][1], hs[k][1]])
        if all(dot(n, x) <= b for n, b in hs):
            verts.add(tuple(x))
    if not verts:
        return None
    x = centroid([list(v) for v in verts])
    if strictly_inside(P, p[-1], A, x) and strictly_inside(P, q[-1], B, x):
        return x
    return None


def det3(m):
    return dot(m[0], cross(m[1], m[2]))


def main():
    name, maxlen = sys.argv[1], int(sys.argv[2])
    r = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    P = Polytope(name)
    U = Unfolder(P, r)
    dist, _ = P.layers(r)
    found = []
    explored = 0
    stack = [(r,)]
    while stack and len(found) < 3:
        path = stack.pop()
        explored += 1
        if len(path) >= 3:
            g_last = cell_geometry(P, path[-1], U.chart3(path))
            for i in range(len(path) - 2):     # skip the parent: it shares a face
                pre = path[:i + 1]
                if certificate(cell_geometry(P, pre[-1], U.chart3(pre)), g_last) is None:
                    found.append((pre, path))
                    break
        if len(path) < maxlen:
            for nb in P.adj[path[-1]]:
                if nb not in path:
                    stack.append(path + (nb,))
    print("%s-cell, root %d, simple paths up to %d cells: explored %d, overlaps found %d"
          % (name, r, maxlen, explored, len(found)))
    for pre, path in found:
        bfs = all(dist[path[i]] == i for i in range(len(path)))
        x = witness(P, U, pre, path)
        print("  cell %d (path %s) overlaps cell %d (path %s)" % (pre[-1], list(pre), path[-1], list(path)))
        print("    path is a BFS path from the root: %s" % bfs)
        print("    common interior point: %s" % ([str(t) for t in x] if x else "not found on the centroid segment"))


if __name__ == "__main__":
    main()
