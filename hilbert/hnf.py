from sympy import pprint, Matrix, floor
from sympy.core.numbers import igcdex


def phi(a, b):
    if b != 0:
        return a - (floor(a/b) * b)
    else:
        return 0


def row_one_gcd(A, i, j, l):
    assert i < A.rows
    assert j < l < A.cols
    if A[i, j] != 0 or A[i, l] != 0:
        u, v, d = igcdex(A[i, j], A[i, l])
        reduced = Matrix([
            [u, -A[i, l]/d],
            [v, A[i, j]/d]
        ])
        s = A.extract(range(A.rows), [j, l]) * (reduced)
        A[:, j] = s[:, 0]
        A[:, l] = s[:, 1]
    return A


def hnf_col(A):
    m, n = A.shape
    A = A.as_mutable().col_join(Matrix.eye(A.cols))
    t = 0
    while t < A.cols:
        r = -1
        s = 0
        while s < A.rows:
            if A[s, r + 1] != 0 or A[s, t] != 0:
                r = r + 1
                if t == r:

                    if A[s, t] < 0:
                        A[:, t] = -1 * A[:, t]
                else:
                    A = row_one_gcd(A, s, r, t)

                for l in range(1, r-1):
                    A[:, l] = A[:, l] - phi(A[s, l], A[s, r]) * A[:, r]
                if t == r:
                    break
            s = s + 1
        t = t + 1

    B = A.extract(range(m), range(n))
    BC = A.extract(range(m, A.rows), range(n))
    return B, BC


def col_one_gcd(A, i, j, l):
    assert i < A.cols
    assert j < l < A.rows
    if A[j, i] != 0 or A[l, i] != 0:
        u, v, d = igcdex(A[j, i], A[l, i])
        print("col_one_gcd")
        print(u, v, d)

        reduced = Matrix([
            [u, v],
            [-A[l, i]/d, A[j, i]/d]
        ])
        s = reduced * A.extract([j, l], range(A.cols))
        # pprint(s)
        A[j, :] = s[0, :]
        A[l, :] = s[1, :]
    return A


def hnf_row(A):
    B, BC = hnf_col(A.T)

    return B.T, BC.T
