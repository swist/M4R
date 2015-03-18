import itertools
from sympy import init_printing, pprint
init_printing(use_unicode=True)


def normal_form(vector_s, set_G, j):
    s_now = vector_s
    while True:
        if any([g == s_now for g in set_G]):
            return s_now - s_now

        normalizable = [g for g in set_G if (g <= s_now)]

        if not normalizable:
            vector_s = s_now
            break
        else:
            vector_g = normalizable[0]

            alpha = s_now.compute_alpha(vector_g)

            s_now = s_now - alpha * vector_g

    return vector_s


def construct_c(G):
    C = []

    # use itertools product instead of zip for better performance
    for f, g in itertools.product(G, G):
        s = f.compute_s_vector(g)

        # only add if well don't have a vector that is equal
        if s and not s.is_zero and all((s != c for c in C)):
            C.append(s)

    return C


def extend_c(f, C, G):

    for g in G:
        vec = f.compute_s_vector(g)

        pprint(vec)
        if vec and all((vec != c for c in C)):
            C.append(vec)

    if all(g != f for g in G):
        G.append(f)


def poittier(F, j):

    G = F
    C = construct_c(G)
    while len(C):
        s = C.pop()

        f = normal_form(s, G, j)

        if f and f.is_zero is False:
            extend_c(f, C, G)

    # return H^+
    return [g for g in G if g[-1] >= 0]
