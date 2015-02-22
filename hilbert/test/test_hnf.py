from sympy import Matrix, gcd, pprint
from hilbert.hnf import row_one_gcd, hnf_col, col_one_gcd, hnf_row


def test_row_one_gcd():
    A = Matrix([[1, 2, 123], [13, 2, 132]])
    B = row_one_gcd(A, 0, 0, 1)

    assert B[0, 0] == gcd(A[0, 0],  A[0, 1])
    assert B[0, 1] == 0


def test_col_one_gcd():
    A = Matrix([[1, 2, 123], [13, 2, 132]]).T

    B = col_one_gcd(A, 0, 1, 2)

    assert B[1, 0] == gcd(A[0, 0],  A[1, 0])
    assert B[2, 0] == 0


def test_hnf_col():
    A = Matrix([[1, 2, 123], [13, 2, 132]])
    B, pivots = hnf_col(A)
    pprint(B)
    assert B.is_lower

    A = Matrix([
        [0, 1, 0],
        [2, -1, 0]
        ])
    B, pivots = hnf_col(A)
    assert B.is_lower
    assert B == Matrix([[1, 0], [0, 2]])


# def test_hnf_row():
#     A = Matrix([
#         [3, 3, 1, 4],
#         [0, 1, 0, 0],
#         [0, 0, 19, 6],
#         [0, 0, 0, 3]
#     ])
#     B, pivots = hnf_row(A)
#     print(B.C)

#     assert B.is_upper
#     assert B == Matrix([
#         [3, 0, 1, 1],
#         [0, 1, 0, 0],
#         [0, 0, 19, 1],
#         [0, 0, 0, 3]
#     ])

#     assert B.is_upper

#     A = Matrix([
#         [0, 1],
#         [2, -1]
#         ])
#     B, pivots = hnf_row(A)
