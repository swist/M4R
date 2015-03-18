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
    A = A.as_mutable()
    i = [0 for i in range(A.cols)]
    t = 0
    while t < A.cols:
        r = -1
        s = 0
        while s < A.rows:
            if A[s, r + 1] != 0 or A[s, t] != 0:
                r = r + 1
                i[r - 1] = s

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
    return A, i


def col_one_gcd(A, i, j, l):
    assert i < A.cols
    assert j < l < A.rows
    if A[j, i] != 0 or A[l, i] != 0:
        u, v, d = igcdex(A[j, i], A[l, i])
        # print(u, v, d)

        reduced = Matrix([
            [u, v],
            [-A[l, i]/d, A[j, i]/d]
        ])
        s = reduced * A.extract([j, l], range(A.cols))
        # pprint(s)
        A[j, :] = s[0, :]
        A[l, :] = s[1, :]
        pprint(A)
    return A


def hnf_row(A):
    # i = [0 for i in range(A.cols)]

    # pprint(A)

    # for t in range(A.rows):
    #     r = -1
    #     for s in range(A.cols):
    #         print("t:%s, s:%s, r:%s" % (t, s, r))
    #         if A[r + 1, s] != 0 or A[t, s] != 0:
    #             r = r + 1
    #             i[r - 1] = s

    #             if t == r:
    #                 if A[t, s] < 0:
    #                     A.row_op(t, lambda x, _: -1*x)
    #             else:
    #                 print("calling with t:%s, s:%s, r:%s" % (t, s, r))
    #                 A = col_one_gcd(A, s, r, t)

    #             for l in range(0, r):
    #                 A[l, :] = A[l, :] - phi(A[r, s], A[l, s]) * A[r, :]
    #             if t == r:
    #                 print(t, r, s)
    #                 print('breaking')
    #                 break
    #         pprint(A)
    # pprint(A)
    # return A, i
    B, BC = hnf_col(A.T)

    return B.T, BC
