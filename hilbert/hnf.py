from sympy import pprint, Matrix, floor
from sympy.core.numbers import igcdex


def phi(a, b):
    if a != 0 and b != 0:
        return a - (floor(a/b) * b)
    else:
        return 1


def row_one_gcd(A, i, j, l):
    assert i < A.rows
    assert j < l < A.cols

    if A[i, j] != 0 or A[i, l] != 0:
        u, v, d = igcdex(A[i, j], A[i, l])
        print(A[i, j], A[i, l])
        reduced = Matrix([
            [u, -A[i, l]/d],
            [v, A[i, j]/d]
        ])

        s = A.extract(range(A.rows), [j, l]) * (reduced)

        A[:, j] = s[:, 0]
        A[:, l] = s[:, 1]
    return A


def hnf_col(A):
    i = [0 for i in range(A.cols)]
    for t in range(A.cols):
        r = -1
        for s in range(A.rows):
            print("t:%s, s:%s, r:%s" % (t, s, r))
            if A[s, r + 1] != 0 or A[s, t] != 0:
                r = r + 1
                i[r - 1] = s
                print("new t:%s, s:%s, r:%s" % (t, s, r))
                pprint(A[s, r])
                pprint(A[s, t])
                if t == r:
                    print("t==r")
                    if A[s, t] < 0:
                        pprint(A[s, t])
                        pprint(A)
                        A[:, t] = -1 * A[:, t]  
                        pprint(A)
                else:
                    A = row_one_gcd(A, s, r, t)

                if r >= 2:
                    for l in range(0, r - 2):
                        A[:, l] = A[:, l] - phi(A[s, r], A[s, l]) * A[:, r]
                        pprint(A[:, l])
                if t == r:
                    break
            pprint(A)
    return A, i


def col_one_gcd(A, i, j, l):

    assert i < A.cols
    assert j < l < A.rows
    pprint(A)
    if A[j, i] != 0 or A[l, i] != 0:
        print(A[j, i])
        print(A[l, i])
        u, v, d = igcdex(A[j, i], A[l, i])
        # print(u, v, d)

        reduced = Matrix([
            [u, v],
            [-A[l, i]/d, A[j, i]/d]
        ])
        pprint(reduced)
        pprint(A.extract([j, l], range(A.cols)))

        s = reduced * A.extract([j, l], range(A.cols))
        # pprint(s)
        A[j, :] = s[0, :]
        A[l, :] = s[1, :]
        pprint(A)
    return A


def hnf_row(A):
    i = [0 for i in range(A.cols)]
    for t in range(A.rows):
        r = -1
        for s in range(A.cols):
            # print("t:%s, s:%s, r:%s" % (t, s, r))
            if A[r + 1, s] != 0 or A[t, s] != 0:
                r = r + 1
                i[r - 1] = s

                if t == r:
                    if A[t, s] < 0:
                        A.col_op(t, lambda x, _: -1*x)
                else:
                    A = col_one_gcd(A, s, r, t)

                for l in range(0, r):
                    A[l, :] = A[l, :] - phi(A[r, s], A[l, s]) * A[r, :]
                if t == r:
                    break
            pprint(A)
    return A, i
