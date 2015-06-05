from sympy import pprint, Matrix
from hilbert.hnf import hnf_row
from hilbert.vector_types import BasisElement
from operator import methodcaller
from orderedset import OrderedSet

def critical_pairs(f, G):
    # Generators automatically use parallel for loops if available
    return {s for s in map(f.s_vector, G) if not s.origin }

def poittier(G, j):
    
    C = reduce(set.union, (critical_pairs(f, G) for f in G))
    pprint(C)
    while len(C):
        s = min(C,key=methodcaller('norm'))
        C.remove(s)
        
        
        if s.irreducible(G):
            G.add(s)
            C |= critical_pairs(s, G)

    return G


def construct_generating_set(A, VectorClass=BasisElement):
    # input - sympy cone in row echelon form
    # pick first generator    
    
    h_1 = VectorClass([A[0, 0]])
    H = {h_1}
    s, n = A.shape
    
    
    for j in range(1, n):
        F = set()
        # project to first j+1 coordinates
        K = A.T[:j+1, :j+1]
        if j < s:
            for h in H:
                F.add(h.lift(K))

            F.add(VectorClass(K[:,-1]))
            F.add(VectorClass(-1*K[:,-1]))
        else:
            for h in H:
                F.add(h.lift(K))
        H = poittier(F, j)
        # pick H^+
        H -= set(filter(lambda x: x[-1] < 0, H))
        pprint(H)

    return H
