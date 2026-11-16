"""known_counterexamples.py -- in the published overlapping ridge unfoldings,
are the two overlapping cells adjacent in P?

In every BFS tree of a regular 4-polytope the overlaps lie across cut ridges
(two cells adjacent in P, not joined in the tree; all_bfs.py).  The known
non-BFS counterexamples:

  seamons24  J. Seamons, Zenodo 2026 (doi:10.5281/zenodo.22004769): the tree
             path [5,14,0,8,10,11,15,3]; cells 5 and 3 overlap.  Facets are
             numbered -e1,+e1,...,-e4,+e4, then the 16 sign vectors (minus
             first); both bit orders are tried.
  dhr24      Devadoss-Harvey-Richter, arXiv:2608.30285, Sec. 2: a chain of 14
             octahedra whose first and last overlap.  Only f1 = (0,2,0,0), M1
             and the unfolded vertices of f2 and f14 are printed, so the chain
             is recovered: every chain from f1 through f2 is unfolded (pruned
             by distance to f14's centroid) and kept iff its 14th cell lands
             exactly on the printed M14 f14.  In their coordinates the 24-cell
             has vertices (+-2,+-2,0,0) permuted, cell centroids (+-2,0,0,0)
             permuted and (+-1,+-1,+-1,+-1), cells adjacent iff C.C' = 2.
  dhr600     their Sec. 4: 8 tetrahedra around one vertex, first and last
             overlapping.  Every 8-cell path in one vertex figure is unfolded
             exactly and the overlapping pairs are listed.
  (120-cell: src/bfsnet/check_dhr_chain_120cell.wls reports it.)

  python3 -I known_counterexamples.py           # all three (about 10 min)
  python3 -I known_counterexamples.py dhr24
"""
import collections, itertools, math, os, sys
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                     # -I drops the script dir

dot = lambda a, b: sum(x * y for x, y in zip(a, b))
sub = lambda a, b: [x - y for x, y in zip(a, b)]
cross = lambda a, b: [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]


def seamons24():
    print("Seamons, 24-cell")
    e = [tuple(s if k == i else 0 for k in range(4)) for i in range(4) for s in (-1, 1)]
    for order in ("most significant bit first", "least significant bit first"):
        sv = []
        for b in range(16):
            bits = [(b >> (3 - k)) & 1 for k in range(4)] if order[0] == "m" else [(b >> k) & 1 for k in range(4)]
            sv.append(tuple(F(1 if t else -1, 2) for t in bits))
        N = e + sv                             # unit facet normals; adjacent iff N.N' = 1/2
        path = [5, 14, 0, 8, 10, 11, 15, 3]
        chain = all(dot(N[a], N[b]) == F(1, 2) for a, b in zip(path, path[1:]))
        print("  %s: path is a chain %s; cell 5 = %s, cell 3 = %s, N5.N3 = %s -> %s"
              % (order, chain, N[5], N[3], dot(N[5], N[3]),
                 "adjacent" if dot(N[5], N[3]) == F(1, 2) else "NOT adjacent"))


def dhr24():
    print("Devadoss-Harvey-Richter, 24-cell")
    V = sorted({p for base in [(2, 2, 0, 0), (2, -2, 0, 0), (-2, -2, 0, 0)] for p in itertools.permutations(base)})
    Cs = list(itertools.product([-1, 1], repeat=4)) + \
        [tuple(s * 2 if i == j else 0 for i in range(4)) for j in range(4) for s in (1, -1)]
    cellV = {c: [v for v in V if dot(c, v) == 4] for c in Cs}
    adj = {c: [d for d in Cs if dot(c, d) == 2] for c in Cs}
    M1 = [[-6, 0, -48, -24], [-48, 0, -6, 24], [24, 216, -24, 42]]
    f1, f2 = (0, 2, 0, 0), (-1, 1, 1, -1)       # f2: the cell holding the triangle they name
    M2f2 = sorted(tuple(4 * x for x in t) for t in
                  [(3, 24, 96), (-24, -3, 96), (12, -12, 87), (-24, 24, 69), (12, 15, 60), (-15, -12, 60)])
    M14f14 = sorted(tuple(4 * x for x in t) for t in
                    [(41, 20, 140), (16, -8, 133), (32, 11, 104), (4, 25, 148), (20, 44, 119), (-5, 16, 112)])

    def child(P, c, d, num):
        """d glued to the placed c: reflect c in the shared triangle; a vertex of
        d is the image of the vertex of c opposite the same triangle vertex."""
        tri = [v for v in cellV[c] if v in cellV[d]]
        a, b, x = (P[v] for v in tri)
        nrm = cross(sub(b, a), sub(x, a))
        nn = dot(nrm, nrm)
        opp = lambda cell, v: next(w for w in cellV[cell] if dot(sub(w, v), sub(w, v)) == 16)
        Q = {v: P[v] for v in tri}
        for u in cellV[c]:
            if u not in tri:
                s = 2 * dot(sub(P[u], a), nrm)
                Q[opp(d, opp(c, u))] = tuple(p - (num(s) / nn) * m for p, m in zip(P[u], nrm))
        return Q

    def place(chain, num):
        P = {v: tuple(num(sum(M1[i][k] * v[k] for k in range(4))) for i in range(3)) for v in cellV[chain[0]]}
        for c, d in zip(chain, chain[1:]):
            P = child(P, c, d, num)
        return P

    cen = lambda P: [sum(x[i] for x in P.values()) / 6 for i in range(3)]
    target = [sum(t[i] for t in M14f14) / 6 for i in range(3)]
    P2 = place([f1, f2], float)
    print("  f2 lands on the printed M2 f2:", sorted(tuple(round(x) for x in p) for p in P2.values()) == M2f2)
    step = math.dist(cen(place([f1], float)), cen(P2))
    found = []

    def dfs(chain, P):
        if math.dist(cen(P), target) > (14 - len(chain)) * step + 1e-6:
            return
        if len(chain) == 14:
            if sorted(tuple(round(x) for x in p) for p in P.values()) == M14f14:
                found.append(chain)
            return
        for d in adj[chain[-1]]:
            if d not in chain:
                dfs(chain + [d], child(P, chain[-1], d, float))
    dfs([f1, f2], P2)
    print("  chains whose 14th cell lands on the printed M14 f14:", len(found))
    for ch in found:
        P = place(ch, F)                            # exact
        exact = sorted(P.values()) == sorted(tuple(F(x) for x in t) for t in M14f14)
        # the overlap they claim: a vertex of f14 strictly inside f1 (f1 = {x : C.x <= 4 on its facets})
        P1 = place(ch[:1], F)
        inside = []
        for w, x in P.items():
            # barycentric test in the 3-flat: x is a convex combination of f1's vertices with
            # x strictly inside iff it lies strictly on the inner side of all 8 face planes of f1
            ok = True
            for c2 in adj[ch[0]]:
                tri = [v for v in cellV[ch[0]] if v in cellV[c2]]
                a, b, y = (P1[v] for v in tri)
                nrm = cross(sub(b, a), sub(y, a))
                other = next(P1[v] for v in cellV[ch[0]] if v not in tri)
                if dot(sub(x, a), nrm) * dot(sub(other, a), nrm) <= 0:
                    ok = False
            if ok:
                inside.append(w)
        f14 = ch[-1]
        print("  ", ch)
        print("    exact placement matches: %s; vertices of f14 strictly inside f1: %d" % (exact, len(inside)))
        print("    f1 = %s, f14 = %s, f1.f14 = %d -> %s, common vertices %d"
              % (ch[0], f14, dot(ch[0], f14), "adjacent" if dot(ch[0], f14) == 2 else "NOT adjacent",
                 len(set(cellV[ch[0]]) & set(cellV[f14]))))


def dhr600():
    print("Devadoss-Harvey-Richter, 600-cell (every 8-cell path in one vertex figure)")
    import certify_bfs_exact as C
    Pt = C.Polytope("600")
    fig = [c for c in range(Pt.nC) if 0 in Pt.cellV[c]]
    fa = {c: [d for d in Pt.adj[c] if d in fig] for c in fig}
    paths = []

    def grow(p):
        if len(p) == 8:
            paths.append(tuple(p))
            return
        for d in fa[p[-1]]:
            if d not in p:
                grow(p + [d])
    for s in fig:
        grow([s])
    paths = sorted({min(p, p[::-1]) for p in paths})
    stat = collections.Counter()
    for p in paths:
        U = C.Unfolder(Pt, p[0])
        g = [C.cell_geometry(Pt, p[i], U.chart3(p[:i + 1])) for i in range(8)]
        ov = [(i, j) for i in range(8) for j in range(i + 2, 8) if C.certificate(g[i], g[j]) is None]
        stat["chains overlapping"] += bool(ov)
        for i, j in ov:
            stat["cells %d and %d, %s" % (i + 1, j + 1, "adjacent in P" if p[j] in Pt.adj[p[i]] else "NOT adjacent")] += 1
    print("  %d paths; %s" % (len(paths), dict(stat)))


if __name__ == "__main__":
    for name in sys.argv[1:] or ["seamons24", "dhr24", "dhr600"]:
        globals()[name]()
