"""
verify_certificate.py
Check a certificate written by certify_bfs_exact.py --export, without
trusting (or importing) any of the code that produced it.

Claim: for the regular 4-polytope P and every root cell r, every BFS
spanning tree rooted at r yields a net whose cells have pairwise
disjoint interiors.

What is checked, all in exact rational arithmetic:

 1. The polytope in the certificate is the one in _4DData/ (the 5-cell,
    whose file uses GoldenRatio, is checked to be conv(e_1..e_5)).
 2. Every placement is the unfolded cell.  The placement of cell c along
    the path (..., p, c) is pinned down by four conditions:
      (a) the vertices of the face shared with p sit where p put them;
      (b) every pairwise distance between vertices of c is preserved;
      (c) every vertex lies in the 3-flat of the root cell;
      (d) the centroids of p and c lie strictly on opposite sides of the
          plane of the shared face within that 3-flat.
    (a)-(c) leave exactly two candidates, mirror images in the face
    plane, and the unfolding rotates c away from p, which is (d).
 3. The placements present are exactly the shortest root paths, and the
    pairs present are exactly the pairs of root paths that agree on the
    parent of every cell they share.  Such a pair occurs together in some
    BFS tree (the other parents are free), and every pair of cells of
    every BFS tree is such a pair -- so these pairs cover every tree.
 4. For each pair a functional u, non-constant on the 3-flat, with
    max_A u.x <= min_B u.x: the two cells lie in opposite closed half
    spaces, so their interiors are disjoint.

Usage: python3 -I verify_certificate.py <cert.json.gz>
"""

import gzip
import itertools
import json
import os
import sys
from fractions import Fraction

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))


def fail(msg):
    print("FAIL:", msg)
    sys.exit(1)


def vec(strs):
    return [Fraction(s) for s in strs]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def minus(a, b):
    return [x - y for x, y in zip(a, b)]


def det(m):
    """Determinant by fraction-exact elimination."""
    m = [list(r) for r in m]
    n, d = len(m), Fraction(1)
    for c in range(n):
        p = next((i for i in range(c, n) if m[i][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            m[c], m[p] = m[p], m[c]
            d = -d
        d *= m[c][c]
        for i in range(c + 1, n):
            f = m[i][c] / m[c][c]
            m[i] = [x - f * y for x, y in zip(m[i], m[c])]
    return d


def rank(rows):
    m = [list(r) for r in rows]
    rk, cols = 0, len(m[0]) if m else 0
    for c in range(cols):
        p = next((i for i in range(rk, len(m)) if m[i][c] != 0), None)
        if p is None:
            continue
        m[rk], m[p] = m[p], m[rk]
        for i in range(len(m)):
            if i != rk and m[i][c] != 0:
                f = m[i][c] / m[rk][c]
                m[i] = [x - f * y for x, y in zip(m[i], m[rk])]
        rk += 1
    return rk


def check_polytope(doc):
    name = doc["polytope"]
    verts = [vec(v) for v in doc["vertices"]]
    faces, cells = doc["faces"], doc["cells"]
    if name == "5":
        if verts != [[Fraction(int(i == j)) for j in range(5)] for i in range(5)]:
            fail("5-cell vertices are not e_1..e_5")
        if sorted(map(sorted, faces)) != [list(t) for t in itertools.combinations(range(5), 3)]:
            fail("5-cell faces are not the 3-subsets")
        cellsets = sorted(sorted({v for f in c for v in faces[f]}) for c in cells)
        if cellsets != [list(t) for t in itertools.combinations(range(5), 4)]:
            fail("5-cell cells are not the 4-subsets")
        for c in cells:
            vs = {v for f in c for v in faces[f]}
            if sorted(c) != [k for k, f in enumerate(faces) if set(f) <= vs]:
                fail("5-cell cell does not list its faces")
    else:
        with open(os.path.join(BASE, "_4DData", "f%s.m" % name)) as fh:
            txt = fh.read().replace("{", "[").replace("}", "]")
        raw = json.loads(txt)
        if verts != [[Fraction(x) for x in v] for v in raw[0]]:
            fail("vertices differ from _4DData")
        if faces != [[i - 1 for i in f] for f in raw[3]] or cells != [[i - 1 for i in c] for c in raw[5]]:
            fail("faces or cells differ from _4DData")
    return verts, faces, cells


def main():
    with gzip.open(sys.argv[1], "rt") as fh:
        doc = json.load(fh)
    verts, faces, cells = check_polytope(doc)
    nC = len(cells)
    cellV = [sorted({v for f in c for v in faces[f]}) for c in cells]
    owners = {}
    for c, fs in enumerate(cells):
        for f in fs:
            owners.setdefault(f, []).append(c)
    if any(len(o) != 2 for o in owners.values()):
        fail("a face does not lie in exactly two cells")
    shared = {}
    adj = [set() for _ in range(nC)]
    for f, (a, b) in owners.items():
        shared[(a, b)] = shared[(b, a)] = f
        adj[a].add(b)
        adj[b].add(a)
    if sorted(rc["root"] for rc in doc["roots"]) != list(range(nC)):
        fail("not every root is present")

    total_trees = total_pairs = 0
    for rc in doc["roots"]:
        r = rc["root"]
        # BFS distances and all shortest root paths
        dist = {r: 0}
        layer = [r]
        while layer:
            nxt = []
            for p in layer:
                for q in adj[p]:
                    if q not in dist:
                        dist[q] = dist[p] + 1
                        nxt.append(q)
            layer = nxt
        paths = {r: [(r,)]}
        for v in sorted(dist, key=dist.get):
            if v != r:
                paths[v] = [pp + (v,) for p in adj[v] if dist[p] == dist[v] - 1 for pp in paths[p]]
        trees = 1
        for v in range(nC):
            if v != r:
                trees *= sum(1 for p in adj[v] if dist[p] == dist[v] - 1)

        # 2. placements
        pl = {tuple(p): {int(k): vec(x) for k, x in m.items()} for p, m in rc["placements"]}
        if set(pl) != {p for ps in paths.values() for p in ps}:
            fail("root %d: placements are not exactly the shortest root paths" % r)
        root = pl[(r,)]
        if root != {v: verts[v] for v in cellV[r]}:
            fail("root %d: root cell is not in place" % r)
        origin = root[cellV[r][0]]
        dirs = [minus(root[v], origin) for v in cellV[r]]
        N = [vec(n) for n in rc["flat_normals"]]
        d = len(origin)
        if rank(dirs) != 3 or len(N) != d - 3 or rank(dirs + N) != d:
            fail("root %d: the flat or its normals are wrong" % r)
        if any(dot(n, x) != 0 for n in N for x in dirs):
            fail("root %d: a normal is not normal to the flat" % r)

        def cent(m):
            return [sum(c) / len(m) for c in zip(*m.values())]

        for path, m in pl.items():
            if len(path) == 1:
                continue
            p, c = path[-2], path[-1]
            par = pl[path[:-1]]
            if sorted(m) != cellV[c]:
                fail("root %d path %s: wrong vertex set" % (r, path))
            fv = faces[shared[(p, c)]]
            if any(m[v] != par[v] for v in fv):                                  # (a)
                fail("root %d path %s: shared face moved" % (r, path))
            for x, y in itertools.combinations(cellV[c], 2):                    # (b)
                if dot(minus(m[x], m[y]), minus(m[x], m[y])) != dot(minus(verts[x], verts[y]), minus(verts[x], verts[y])):
                    fail("root %d path %s: not an isometry" % (r, path))
            if any(dot(n, minus(x, origin)) != 0 for n in N for x in m.values()):  # (c)
                fail("root %d path %s: left the root flat" % (r, path))
            F = [m[v] for v in fv]
            F0 = F[0]
            E = [minus(x, F0) for x in F[1:]]
            E = next([a, b] for a, b in itertools.combinations(E, 2) if rank([a, b]) == 2)

            def side(x):
                return det(E + [minus(x, F0)] + N)

            if not side(cent(par)) * side(cent(m)) < 0:                          # (d)
                fail("root %d path %s: not on the far side of the shared face" % (r, path))

        # 3. the pairs present are exactly the compatible pairs of root paths
        def compatible(P, Q):
            mp = {P[i]: P[i - 1] for i in range(1, len(P))}
            return all(mp.get(Q[i], Q[i - 1]) == Q[i - 1] for i in range(1, len(Q)))

        need = {(p, q) for u, v in itertools.combinations(range(nC), 2)
                for p in paths[u] for q in paths[v] if compatible(p, q)}
        got = {}
        for p, q, u in rc["pairs"]:
            got[(tuple(p), tuple(q))] = vec(u)
        if set(got) != need:
            fail("root %d: pairs present are not exactly the compatible pairs (%d vs %d)"
                 % (r, len(got), len(need)))

        # 4. separation
        for (p, q), u in got.items():
            if all(dot(u, x) == 0 for x in dirs):
                fail("root %d: functional constant on the flat" % r)
            if max(dot(u, x) for x in pl[p].values()) > min(dot(u, x) for x in pl[q].values()):
                fail("root %d pair %s %s: not separated" % (r, p, q))
        total_trees += trees
        total_pairs += len(got)

    print("OK  %s-cell: %d roots, %d BFS trees, %d placements checked, %d pairs separated"
          % (doc["polytope"], nC, total_trees, sum(len(rc["placements"]) for rc in doc["roots"]), total_pairs))


if __name__ == "__main__":
    main()
