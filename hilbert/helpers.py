import itertools
from sympy import init_printing, pprint
init_printing(use_unicode=True)


def normal_form(vector_s, set_G, j):
    print("computing normal form from:")
    print("s:")
    pprint(vector_s)
    print("G:")
    pprint(set_G)
    s_now = vector_s
    while True:
        if any([g == s_now for g in set_G]):
            return s_now - s_now

        normalizable = [g for g in set_G if (g <= s_now)]
        print("normalizable")
        pprint(normalizable)
        pprint(s_now)
        if not normalizable:
            vector_s = s_now
            break
        else:
            vector_g = normalizable[0]
            print("reducing s:")
            alpha = s_now.compute_alpha(vector_g)
            print("alpha")
            print(alpha)
            s_now = s_now - alpha * vector_g

    return vector_s


def construct_c(G):
    C = []
    print("constructing C:\n")
    # use itertools product instead of zip for better performance
    for f, g in itertools.product(G, G):
        s = f.compute_s_vector(g)

        # only add if well don't have a vector that is equal
        if s and not s.is_zero and all((s != c for c in C)):
            C.append(s)
    print("C is:\n")
    pprint(C)
    return C


def extend_c(f, C, G):
    print("f not zero, extending C and G:\n")
    for g in G:
        vec = f.compute_s_vector(g)
        print("computed s-vector:")
        pprint(vec)
        if vec and all((vec != c for c in C)):
            C.append(vec)

    if all(g != f for g in G):
        G.append(f)


def poittier(F, j):
    print("reduction algorithm:")
    G = F
    C = construct_c(G)
    print("reduction loop:")
    while len(C):

        print("C is:")
        pprint(C)
        print("G is:")
        pprint(G)
        s = C.pop()
        print("s is:")
        pprint(s)

        f = normal_form(s, G, j)

        print("f is:")
        pprint(f)
        if f and f.is_zero is False:
            extend_c(f, C, G)
        print("end of iteration \n")
    # return H^+
    return [g for g in G if g[-1] >= 0]


def ref(matrix):
    rows, cols = matrix.shape
    M = matrix.as_mutable()
    m = min(cols, rows)
    i, stopper = 0, 0
    BC = eye(max(rows, cols))
    while i < m and stopper < cols:
        j = next((k for k in xrange(i, rows) if M[k, i] != 0), None)
        if j is not None:
            for k in xrange(j+1, rows):
                if M[k, i] != 0:
                    a = M[j, i]
                    b = M[k, i]
                    L = gcdex(a, b)
                    row_j = L[0]*M[j, :] + L[1]*M[k, :]
                    row_k = -b/L[2] * M[j, :] + a/L[2] * M[k, :]
                    M[j, :] = row_j
                    M[k, :] = row_k
            if M[i, i] < 0:
                M.row_op(i, lambda x, _: -1*x)
        else:
            M.col_swap(i, cols-1)
            i = i - 1
        i = i + 1
        stopper = stopper + 1

    return matrix._new(M), BC
