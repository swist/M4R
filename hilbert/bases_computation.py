from sympy import pprint, Matrix
from hilbert.helpers import poittier
from hilbert.hnf import hnf_row
from hilbert.vector_types import BasisElement


def construct_generating_set(A, VectorClass=BasisElement):
    # input - sympy cone in row echelon form
    # pick first generator    
    # A = matrix
    h_1 = VectorClass(Matrix([A[0, 0]]))
    H = [h_1]
    s, n = A.shape
    print("A")
    pprint(A)
    # pprint(BC)

    for j in range(1, n):
        
        # pick H^+
        H = [h for h in H if h[-1] >= 0]

        F = []
        # project to first j+1 coordinates
        K_j_1 = A[:j+1, :j + 1]

        if j < s:
            for h in H:
                F.append(h.lift_multiple_choice(K_j_1))

            F.append(VectorClass(A[j, :j+1]))
            F.append(VectorClass(-1*A[j, :j+1]))
        else:
            # print "j>=s \n"
            for h in H:
                F.append(h.lift_single_choice(K_j_1))

        H = poittier(F, j)

    return H
