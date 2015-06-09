from sympy import Matrix, gcd, pprint, randMatrix, eye
from hilbert.hnf import row_one_gcd, hnf_col, hnf_row


def test_row_one_gcd():
    A = Matrix([[1, 2, 123], [13, 2, 132]])
    B = row_one_gcd(A, 0, 0, 1)

    assert B[0, 0] == gcd(A[0, 0],  A[0, 1])
    assert B[0, 1] == 0

def test_hnf_equiv():
    A = Matrix([[1, 2, 123], [13, 2, 132]])
    A, BA = hnf_col(A.T)
    B, BC = hnf_row(A.T)
    pprint(A)
    pprint(B.T)
    assert A == B.T


def test_hnf_col():
    A = Matrix([[1, 2, 123], [13, 2, 132]])
    B, pivots = hnf_col(A)
    pprint(B)
    assert B.is_lower

    A = Matrix([
        [0, 1],
        [2, -1]
        ])
    B, pivots = hnf_col(A)
    assert B.is_lower
    assert B == Matrix([
        [1, 0],
        [-1, 2]
    ])


def test_hnf_row():
    A = Matrix([
        [2, 0],
        [0, 1]
    ])
    B, pivots = hnf_row(A)

    assert B.is_upper
    assert B == A

    A1 = Matrix([
        [2, -1],
        [0, 1]
    ])
    B, pivots = hnf_row(A1)
    assert B.is_upper
    assert B == A1

    A = Matrix([[1, 2, 123], [13, 2, 132]])
    B, pivots = hnf_row(A)
    pprint(B)
    assert B == Matrix([[1, 2, 123], [0, 24, 1467]])

    A = Matrix([[2, 0, -1], [1, 3, 2], [1, 4, 8]])
    B, pivots = hnf_row(A)

    assert B == Matrix([[1, 4, 8], [0, 1, 6], [0, 0, 31]])
