"""tree_observations_120.py -- three observations on the BFS trees of the 120-cell.

A tree's overlaps are its straight 7-cycles (claim 4 of README.md), so they
are counted here combinatorially, without unfolding.

  1. Direction.  Build trees layer by layer, each cell preferring a parent
     through which its grandparent -> parent -> cell step is bent (or, the
     other way, straight).  Preferring bent steps gives valid trees.
  2. Rings.  A straight step x -> m -> v runs along one of the 10-cell
     great-circle rings of the 120-cell; the ring either passes through the
     root (radial) or misses it by one cell.  The straight branch of every
     overlapping 7-cycle misses the root; its partner branch is bent and
     its ring through the root.
  3. Shape.  Valid and invalid trees can be isomorphic as rooted trees, so
     the shape of a tree does not decide validity; its embedding in G(P) does.

  python3 -I tree_observations_120.py        # about 1 min
"""
import collections, gzip, json, math, os, random, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                     # -I drops the script dir
from certify_bfs_exact import Polytope
DATA = os.path.join(HERE, "..", "..", "data", "certify")
P = Polytope("120")
V = [set(s) for s in P.cellV]
ADJ = [set(a) for a in P.adj]
LAY = {}


def layers(r):
    if r not in LAY:
        LAY[r] = P.layers(r)
    return LAY[r]


def straight(x, m, q):
    """x -> m -> q passes straight through m: faces m&x and m&q are opposite."""
    F1, F2 = V[x] & V[m], V[q] & V[m]
    return not any((s in F1 and t in F2) or (t in F1 and s in F2) for s, t in P.cellE[m])


def cycles7(r, parent):
    """The same-layer 7-cycles around a vertex, as (chord branch, other branch)."""
    dist = layers(r)[0]

    def path(v):
        o = [v]
        while o[-1] != r:
            o.append(parent[o[-1]])
        return o[::-1]
    out = []
    for a in range(P.nC):
        for b in P.adj[a]:
            if b <= a or dist[a] != dist[b]:
                continue
            pa, pb = parent[a], parent[b]
            if pa == pb or pb not in ADJ[pa] or b not in ADJ[pa] or a not in ADJ[pb]:
                continue
            A, B = path(a), path(b)
            if len(A) < 4 or A[-4] != B[-4] or A[-3] == B[-3]:
                continue
            x, y = A[-3], B[-3]
            out.append(((x, pa, a), (y, pb, b)) if pa in ADJ[y] else ((y, pb, b), (x, pa, a)))
    return out


def overlaps(r, parent):
    return sum(straight(*side) for side, _ in cycles7(r, parent))


# rings: a pair of adjacent cells spans a 2-plane through the origin; the ring
# is the set of cells whose centroids lie in it
C = []
for c in P.cellCent:
    w = [float(x) for x in c]
    s = math.sqrt(sum(t * t for t in w))
    C.append([t / s for t in w])
dot = lambda a, b: sum(x * y for x, y in zip(a, b))
RINGS = {}


def ring(x, m):
    key = (min(x, m), max(x, m))
    if key not in RINGS:
        e1 = C[x]
        t = [p - dot(C[m], e1) * q for p, q in zip(C[m], e1)]
        s = math.sqrt(dot(t, t))
        e2 = [p / s for p in t]
        RINGS[key] = frozenset(c for c in range(P.nC) if 1 - dot(C[c], e1) ** 2 - dot(C[c], e2) ** 2 < 1e-9)
    return RINGS[key]


def ring_miss(r, x, m):
    """Graph distance from the root to the ring through x and m (0: radial)."""
    d = layers(r)[0]
    return min(d[c] for c in ring(x, m))


def build(r, prefer_straight, rng):
    dist, opts = layers(r)
    parent = [-1] * P.nC
    for v in sorted(range(P.nC), key=lambda c: dist[c]):
        if v == r:
            continue
        o = opts[v]
        if dist[v] >= 2:
            o = [p for p in o if straight(parent[p], p, v) == prefer_straight] or o
        parent[v] = rng.choice(o)
    return parent


def queue_tree(r):
    return [-1 if p is None else p for p in P.queue_tree(r)]


def shape(parent, r):
    ch = collections.defaultdict(list)
    for v, p in enumerate(parent):
        if p >= 0:
            ch[p].append(v)

    def f(v):
        return "(" + "".join(sorted(f(c) for c in ch[v])) + ")"
    return f(r)


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    load = lambda f: json.load(gzip.open(os.path.join(DATA, f), "rt"))
    fams = {"T_r": [(r, queue_tree(r)) for r in range(P.nC)]}
    rng = random.Random(1)
    print("1. direction preference (20 roots x 10 trees each)")
    for prefer, lab in [(False, "prefer bent"), (True, "prefer straight")]:
        trees = [(r, build(r, prefer, rng)) for r in range(0, 120, 6) for _ in range(10)]
        print("   %-16s valid %3d of %d" % (lab, sum(overlaps(r, p) == 0 for r, p in trees), len(trees)))
        if not prefer:
            fams["prefer bent"] = trees
    fams["random"] = [(r, p) for r, ts in enumerate(load("random_trees_120.json.gz")) for p in ts]
    fams["sector"] = [(0, p) for p in load("sector_trees_120.json.gz")[0]]

    print("2. rings of straight steps: distance from the root to the ring (0 = radial)")
    for lab, trees in fams.items():
        steps, chord, other = collections.Counter(), collections.Counter(), collections.Counter()
        for r, p in trees:
            for v in range(P.nC):
                if v != r and p[v] != r and straight(p[p[v]], p[v], v):
                    steps[ring_miss(r, p[p[v]], p[v])] += 1
            for side, oth in cycles7(r, p):
                if straight(*side):
                    chord[ring_miss(r, side[0], side[1])] += 1
                    other[("straight" if straight(*oth) else "bent", ring_miss(r, oth[0], oth[1]))] += 1
        print("   %-12s %5d trees  all straight steps %s  overlapping chord branches %s  their partners %s"
              % (lab, len(trees), dict(sorted(steps.items())), dict(chord), dict(other)))

    print("3. rooted-tree shapes")
    by = collections.defaultdict(lambda: [0, 0])
    for lab, trees in fams.items():
        for r, p in trees:
            by[shape(p, r)][overlaps(r, p) > 0] += 1
    mixed = [k for k, (g, b) in by.items() if g and b]
    print("   %d trees, %d shapes, %d shapes with both valid and invalid trees (%d valid, %d invalid)"
          % (sum(g + b for g, b in by.values()), len(by), len(mixed),
             sum(by[k][0] for k in mixed), sum(by[k][1] for k in mixed)))
