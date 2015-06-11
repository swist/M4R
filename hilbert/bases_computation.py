from sympy import pprint, Matrix
from hilbert.hnf import hnf_row
from hilbert.vector_types import BasisElement
from operator import methodcaller



def preimage(matrix, vectors, VectorClass):
    nullspaces = [(-1 * v).row_join(matrix).nullspace() for v in vectors]
    
    nullspaces = [to_postive_halfspace(ns[0]) for ns in nullspaces]
    vectors = [VectorClass(ns[1:matrix.cols+1]) for ns in nullspaces]
    
    return set(map(convert_to_ZZ, vectors))

def to_postive_halfspace(vector):
    # make sure 1 in the first position
    if vector[0] < 0:
        vector = vector * abs(vector[0])/vector[0]
    return vector

def convert_to_ZZ(vector):
    
    for i in range(len(vector)):
        if(vector[i].is_Rational):
            vector = vector[i].q * vector
    
    return vector

def critical_pairs(s, G):
    return {v for v in map(s.s_vector, G) if not v.origin }

def is_irreducible(s, G):
    return not any(g <= s for g in G)

def construct_generating_set(A, VectorClass=BasisElement):
    H = {VectorClass([A[0, 0]])}
    s, n = A.shape

    for j in range(1, n):
        K = A.T[:j+1, :j+1]
        F = set()
        if j < s:
            F.add(VectorClass(K[:,-1]))
            F.add(VectorClass(-1*K[:,-1]))
        for h in H:
            F.add(h.lift(K))
                
        H = cpc(F)
        # pick H^+
        H -= set(filter(lambda x: x[-1] < 0, H))

    return H

def cpc(G):
    C = reduce(set.union, (critical_pairs(f, G) for f in G))

    while len(C):

        s = min(C,key=methodcaller('norm'))
        C.remove(s)

        if is_irreducible(s, G):
            G.add(s)
            C |= critical_pairs(s, G)

    return G
