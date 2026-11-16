"""path_classes.py -- the root-path pairs of a regular 4-polytope up to symmetry,
and the holonomy of every cut ridge.

Two cells u, v of a BFS tree are placed in the net by their tree path
u ... L ... v alone, where L is the cell at which their root paths part.  So
all_bfs.py's pairs of compatible root paths (from root 0) fall into classes:
tree paths that a symmetry of P (or reversal) carries onto each other give
congruent placements and the same verdict.  This script

  1. builds the symmetry group as permutations of the cells (orthogonal maps
     of the cell centroids; |G| = 384 / 1152 / 14400 for the 16 / 24 / 120-
     and 600-cell),
  2. sorts the pairs into classes and checks that every class is either all
     overlapping or all disjoint (against data/certify/all_bfs_N_overlaps.json.gz;
     for the 16- and 24-cell claim 1 of README.md says no pair overlaps),
  3. for every class of a cut ridge -- u, v adjacent in P but not joined by a
     tree edge -- computes the holonomy of its fundamental cycle: the rigid
     motion of R^3 that takes v from where it would sit glued to u (u's face
     reflected) to where the tree puts it.  Reported: the rotation angle, the
     distance between the centroids of u and v relative to that when they
     touch, and v's displacement along u's outward face normal (> 0: away
     from u) and across it, both relative to the same distance.

Floating point is used for the holonomy only; the verdicts come from the
exact runs of all_bfs.py.

  python3 -I path_classes.py 120            # summary by cycle length (13 s + group)
  python3 -I path_classes.py 120 --all      # one line per cut-ridge class
"""
import collections, gzip, itertools, json, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                     # -I drops the script dir
import certify_bfs_exact as C

NAME = sys.argv[1]
Pt = C.Polytope(NAME)
n, r = Pt.nC, 0
dist, opts = Pt.layers(r)
paths = C.root_paths(r, dist, opts)
ADJ = [set(a) for a in Pt.adj]


# ------------------------------------------------------- symmetry group

def dot(a, b): return sum(x * y for x, y in zip(a, b))
def sub(a, b): return [x - y for x, y in zip(a, b)]
def cross(a, b): return [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
def unit(a): s = math.sqrt(dot(a, a)); return [x / s for x in a]
def mean(ps): return [sum(p[i] for p in ps) / len(ps) for i in range(len(ps[0]))]


def solve(M, B):
    """M X = B for square M (Gaussian elimination, partial pivoting)."""
    k = len(M)
    A = [list(M[i]) + list(B[i]) for i in range(k)]
    for c in range(k):
        p = max(range(c, k), key=lambda i: abs(A[i][c]))
        A[c], A[p] = A[p], A[c]
        if abs(A[c][c]) < 1e-12:
            return None
        A[c] = [x / A[c][c] for x in A[c]]
        for i in range(k):
            if i != c:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [row[k:] for row in A]


def symmetry_group():
    """Cell permutations induced by the orthogonal maps fixing the centroids.
    An orthogonal map is fixed by the images of cell 0 and three neighbours
    spanning R^4 with it; each candidate is kept iff it permutes the cells."""
    X = [[float(x) for x in c] for c in Pt.cellCent]
    key = lambda p: tuple(round(t * 1e6) for t in p)
    index = {key(p): i for i, p in enumerate(X)}
    d2 = lambda a, b: dot(sub(a, b), sub(a, b))
    for t3 in itertools.combinations(Pt.adj[0], 3):
        src = [X[0]] + [X[t] for t in t3]
        if solve(src, [[float(i == j) for j in range(4)] for i in range(4)]) is not None:
            break
    G = []
    for t in range(n):
        for m3 in itertools.permutations(Pt.adj[t], 3):
            if any(abs(d2(X[m3[i]], X[m3[j]]) - d2(X[t3[i]], X[t3[j]])) > 1e-6
                   for i, j in [(0, 1), (0, 2), (1, 2)]):
                continue
            A = solve(src, [X[t]] + [X[m] for m in m3])     # rows: x^T A = image^T
            perm = []
            for p in X:
                j = index.get(key([sum(p[k] * A[k][i] for k in range(4)) for i in range(4)]))
                if j is None:
                    break
                perm.append(j)
            if len(perm) == n:
                G.append(perm)
    return G


G = symmetry_group()
to0 = collections.defaultdict(list)
for g in G:
    to0[g.index(0)].append(g)


def canon(S):
    """Lexicographically least image of the cell sequence S under G and reversal.
    The least image starts with cell 0 (G is cell-transitive), and the
    candidates are pruned one position at a time."""
    best = None
    for T in (S, S[::-1]):
        cand = to0[T[0]]
        for i in range(1, len(T)):
            vals = [g[T[i]] for g in cand]
            m = min(vals)
            cand = [g for g, v in zip(cand, vals) if v == m]
        img = tuple(cand[0][x] for x in T)
        if best is None or img < best:
            best = img
    return best


# ------------------------------------------------------------- classes

ovf = os.path.join(HERE, "..", "..", "data", "certify", "all_bfs_%s_overlaps.json.gz" % NAME)
overlapping = ({(tuple(p), tuple(q)) for p, q in json.load(gzip.open(ovf, "rt"))}
               if os.path.exists(ovf) else set())
classes = {}                                  # canon -> [pairs, verdicts, representative]
for u, v in itertools.combinations(range(n), 2):
    for p in paths[u]:
        for q in paths[v]:
            if not C.compatible(p, q):
                continue
            k = 0
            while k < min(len(p), len(q)) and p[k] == q[k]:
                k += 1
            S = tuple(p[k - 1:][::-1]) + tuple(q[k:])      # the tree path u ... L ... v
            e = classes.setdefault(canon(S), [0, set(), (p, q)])
            e[0] += 1
            e[1].add((p, q) in overlapping)
assert all(len(e[1]) == 1 for e in classes.values()), "a class with both verdicts"


# ------------------------------------------------------------ holonomy

U = C.Unfolder(Pt, r)
B = []                                         # orthonormal basis of the root 3-flat
for b in U.basis:
    w = [float(x) for x in b]
    for e in B:
        w = sub(w, [dot(w, e) * t for t in e])
    B.append(unit(w))
O = [float(x) for x in U.origin]
loc = lambda x: [dot([float(t) - o for t, o in zip(x, O)], e) for e in B]


def holonomy(p, q):
    u, v = p[-1], q[-1]
    glued, tree, cu = U.place(p + (v,)), U.place(q), U.place(p)
    W = sorted(glued)
    Pg, Pt_ = [loc(glued[w]) for w in W], [loc(tree[w]) for w in W]
    for tri in itertools.combinations(range(1, len(W)), 3):     # linear part
        Dg = [sub(Pg[j], Pg[0]) for j in tri]
        if abs(dot(Dg[0], cross(Dg[1], Dg[2]))) > 1e-6:
            break
    Dt = [sub(Pt_[j], Pt_[0]) for j in tri]
    Mt = solve(Dg, Dt)                          # rows: d_glued^T Mt = d_tree^T
    angle = math.degrees(math.acos(max(-1, min(1, (Mt[0][0] + Mt[1][1] + Mt[2][2] - 1) / 2))))
    cg, ct, ca = mean(Pg), mean(Pt_), mean([loc(x) for x in cu.values()])
    fc = mean([loc(cu[w]) for w in Pt.faces[Pt.shared[(u, v)]]])
    touch = math.dist(ca, cg)
    move = sub(ct, cg)
    normal = dot(move, unit(sub(fc, ca)))
    across = math.sqrt(max(0, dot(move, move) - normal ** 2))
    return angle, math.dist(ca, ct) / touch, normal / touch, across / touch


rows = []
kinds = collections.Counter()
for c, (cnt, verdict, (p, q)) in classes.items():
    u, v = p[-1], q[-1]
    k = 0
    while k < min(len(p), len(q)) and p[k] == q[k]:
        k += 1
    if v not in ADJ[u]:
        kinds["not adjacent in P"] += 1
        continue
    if k == min(len(p), len(q)):
        kinds["tree edge"] += 1
        continue
    kinds["cut ridge"] += 1
    length = (len(p) - k) + (len(q) - k) + 1
    rows.append(("same" if dist[u] == dist[v] else "next", length) + holonomy(p, q)
                + (True in verdict, cnt))

print("%s-cell: |G| = %d, %d root-path pairs from root 0, %d classes: %s"
      % (NAME, len(G), sum(e[0] for e in classes.values()), len(classes),
         ", ".join("%s %d" % kv for kv in sorted(kinds.items()))))
print("overlapping classes: %d" % sum(True in e[1] for e in classes.values()))
print()
print("cut ridges by fundamental cycle (layers of u, v; cycle length):")
print("layer len classes   rotation (deg)   min dist/touch   min normal   across at min   overlap")
for key in sorted({row[:2] for row in rows}):
    R = [row for row in rows if row[:2] == key]
    m = min(R, key=lambda row: row[3])
    print("%-5s %3d %6d   %7.2f..%-7.2f   %8.3f        %+8.3f      %7.3f        %s"
          % (key + (len(R), min(row[2] for row in R), max(row[2] for row in R), m[3],
                    min(row[4] for row in R), m[5], any(row[6] for row in R))))
if "--all" in sys.argv:
    print()
    print("layer len  rotation  dist/touch  normal   across  overlap  pairs")
    for row in sorted(rows):
        print("%-5s %3d  %7.2f   %7.3f   %+7.3f  %7.3f  %-7s %6d" % row)
