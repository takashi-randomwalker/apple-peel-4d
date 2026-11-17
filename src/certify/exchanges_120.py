"""exchanges_120.py -- exchanges in the BFS trees of the 120-cell.

An exchange is a pair of adjacent cells a, b of one layer with different
parents p_a != p_b, where p_b is adjacent to a and p_a to b: each could have
taken the other's parent.  Every straight 7-cycle contains one, and by
the paper's lemma "What a traversal can reach" no tree produced by a
breadth-first traversal does (it would need p_a before p_b and p_b before
p_a in the queue order) -- this is the corollary "Traversal trees".  This script counts them in T_r from
every root and in the random trees of data/certify/random_trees_120.json.gz.

  python3 -I exchanges_120.py
"""
import gzip, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                     # -I drops the script dir
from certify_bfs_exact import Polytope

P = Polytope("120")
ADJ = [set(a) for a in P.adj]


def exchanges(r, parent):
    d = P.layers(r)[0]
    return sum(1 for a in range(P.nC) for b in P.adj[a]
               if b > a and a != r and d[a] == d[b] and parent[a] != parent[b]
               and parent[b] in ADJ[a] and parent[a] in ADJ[b])


tq = [exchanges(r, P.queue_tree(r)) for r in range(P.nC)]
print("T_r, 120 roots: %d exchanges" % sum(tq))
rnd = json.load(gzip.open(os.path.join(HERE, "..", "..", "data", "certify",
                                       "random_trees_120.json.gz"), "rt"))
c = [exchanges(r, p) for r, ts in enumerate(rnd) for p in ts]
print("random trees: %d, with an exchange %d, exchanges in all %d (min %d, max %d per tree)"
      % (len(c), sum(x > 0 for x in c), sum(c), min(c), max(c)))
