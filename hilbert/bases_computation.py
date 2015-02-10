import numpy as np
import math
from sympy import *
from .helpers import ref, poittier
from .vector_types import BasisElement, ExtremeRay


def construct_hilbert_basis(matrix, VectorClass=BasisElement):
    # input - sympy matrix in row echelon form
    # pick first generator
    pprint(VectorClass)
    A, BC = ref(matrix)
    # A = matrix
    h_1 = VectorClass(Matrix([A[0, 0]]))
    H = [h_1]
    s, n = A.shape
    print("A, BC")
    pprint(A)
    # pprint(BC)

    for j in xrange(1, n):
        # print 'H'
        # pprint(H)
        F = []
        # project to first j+1 coordinates
        K_j_1 = A[:j+1, j]
        print("H_%d to be lifted:" % j)
        pprint(H)

        if j < s:
            for h in H:
                F.append(h.lift_multiple_choice(K_j_1))

            F.append(VectorClass(A[j, :j+1]))
            F.append(VectorClass(-1*A[j, :j+1]))
        else:
            # print "j>=s \n"
            for h in H:
                F.append(h.lift_single_choice(K_j_1))

        print("F after lifting:")
        pprint(F)

        H = poittier(F, j)
        print("H_%d:" % (j+1))
        pprint(H)
    return H
