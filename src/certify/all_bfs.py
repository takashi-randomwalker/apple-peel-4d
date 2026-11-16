"""all_bfs.py -- every BFS tree of the 120- or 600-cell, exactly.

A cell's placement in the net depends only on its path to the root, so a
BFS tree gives a net iff no two of its root paths overlap.  All pairs of
root paths from root 0 that occur together in some BFS tree (they give
every shared cell the same parent) are decided exactly over Q(sqrt5) with
certify_bfs_exact.certificate.  Both polytopes are cell-transitive, so
root 0 stands for every root.

  python3 -I all_bfs.py 120            # certify all pairs, write
                                       # data/certify/all_bfs_120_overlaps.json.gz
  python3 -I all_bfs.py 120 classify   # sort the pairs by cycle type (seconds)

120-cell (2026-11-16): 901 root paths, 388,464 pairs, 53 min with 7
processes; 240 overlap, all of one type -- see classify.
600-cell: 1,889 root paths, 1,757,416 pairs."""
import collections, gzip, itertools, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                     # -I drops the script dir
import certify_bfs_exact as C
NAME = sys.argv[1] if len(sys.argv) > 1 else "120"
OUT = os.path.join(HERE, "..", "..", "data", "certify", "all_bfs_%s_overlaps.json.gz" % NAME)
r = 0
Pt = C.Polytope(NAME)
dist, opts = Pt.layers(r)
paths = C.root_paths(r, dist, opts)
ADJ = [set(a) for a in Pt.adj]
V = [set(s) for s in Pt.cellV]


def work(u):
    U = C.Unfolder(Pt, r)
    geo = {}

    def g(p):
        if p not in geo:
            geo[p] = C.cell_geometry(Pt, p[-1], U.chart3(p))
        return geo[p]
    n, ov = 0, []
    for v in range(u + 1, Pt.nC):
        for p in paths[u]:
            for q in paths[v]:
                if C.compatible(p, q):
                    n += 1
                    if C.certificate(g(p), g(q)) is None:
                        ov.append((p, q))
    return u, n, ov


def straight(x, m, q):
    """x -> m -> q passes straight through m: faces m&x and m&q are opposite."""
    F1, F2 = V[x] & V[m], V[q] & V[m]
    return not any((s in F1 and t in F2) or (t in F1 and s in F2) for s, t in Pt.cellE[m])


def classify():
    ov = {(tuple(p), tuple(q)) for p, q in json.load(gzip.open(OUT, "rt"))}
    cnt = collections.Counter()
    seen = 0
    for u, v in itertools.combinations(range(Pt.nC), 2):
        for p in paths[u]:
            for q in paths[v]:
                if not C.compatible(p, q):
                    continue
                k = 0
                while k < min(len(p), len(q)) and p[k] == q[k]:
                    k += 1
                length = (len(p) - k) + (len(q) - k) + 1
                hit = (p, q) in ov
                seen += hit
                if v not in ADJ[u]:
                    cnt[("not adjacent in P", hit)] += 1
                    continue
                if k == min(len(p), len(q)):                     # one path extends the other
                    cnt[("adjacent, joined by a tree edge", hit)] += 1
                    continue
                if dist[u] != dist[v] or length != 7:
                    cnt[("adjacent, %s layer, cycle %d" % ("same" if dist[u] == dist[v] else "next", length), hit)] += 1
                    continue
                pa, pb, x, y = p[-2], q[-2], p[-3], q[-3]
                if not (pb in ADJ[pa] and v in ADJ[pa] and u in ADJ[pb]):
                    cnt[("same layer, cycle 7, {a,b,pa,pb} not a clique", hit)] += 1
                    continue
                if (pa in ADJ[y]) == (pb in ADJ[x]):            # the 120-cell always has exactly one chord
                    cnt[("{a,b,pa,pb} a clique, %s chords" % ("two" if pa in ADJ[y] else "no"), hit)] += 1
                    continue
                side, other = ((x, pa, u), (y, pb, v)) if pa in ADJ[y] else ((y, pb, v), (x, pa, u))
                cnt[("{a,b,pa,pb} a clique: chord side %s, other side %s"
                     % ("straight" if straight(*side) else "bent",
                        "straight" if straight(*other) else "bent"), hit)] += 1
    assert seen == len(ov)
    for (k, hit), n in sorted(cnt.items()):
        print("  %-62s %-8s %7d" % (k, "OVERLAP" if hit else "ok", n))


if __name__ == "__main__":
    if sys.argv[2:] == ["classify"]:
        classify()
        sys.exit()
    from concurrent.futures import ProcessPoolExecutor
    t0, tot, allov = time.time(), 0, []
    with ProcessPoolExecutor(7) as ex:
        for u, n, ov in ex.map(work, sorted(range(Pt.nC), key=lambda u: -len(paths[u]))):
            tot += n
            allov += ov
            print("cell %3d pairs %7d overlaps %5d  (%.0fs)" % (u, n, len(ov), time.time() - t0), flush=True)
    print("TOTAL pairs", tot, "overlapping", len(allov), "(%.0fs)" % (time.time() - t0))
    with gzip.open(OUT, "wt") as fh:
        json.dump([[list(p), list(q)] for p, q in allov], fh)
