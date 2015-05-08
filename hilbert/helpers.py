import itertools
from operator import methodcaller
from sympy import init_printing, pprint
init_printing(use_unicode=True)


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
        if vec and all((vec != c for c in C)):
            C.append(vec)

    if all(g != f for g in G):
        G.append(f)


def normal_form(vector_s, set_G, j):
    s_now = vector_s
    while True:
        if any([g == s_now for g in set_G]):
            return s_now - s_now

        normalizable = [g for g in set_G if (g <= s_now)]
        print('normalizing')
        pprint(normalizable)
        if not normalizable:
            vector_s = s_now
            break
        else:
            vector_g = normalizable[0]
            alpha = s_now.compute_alpha(vector_g)
            s_now = s_now - alpha * vector_g

    return vector_s


def poittier2(F, j):
    G = F
    C = construct_c(G)

    while len(C):
        print("G, n.els:%d" %len(G)) 
        pprint(G)
        print("\nC, n.els:%d" %len(C)) 
        pprint(C)
        s = C.pop()
        pprint([g <= s for g in G])
        f = normal_form(s, G, j)

        if f and f.is_zero is False:
            extend_c(f, C, G)

    # return H^+
    return [g for g in G if g[-1] >= 0]


def poittier(G, j):
    
    C = construct_c(G)

    while len(C):
        C = sorted(C,key=methodcaller('norm'), reverse=True)
        print("G, n.els:%d" %len(G)) 
        pprint(G)
        print("\nC, n.els:%d" %len(C)) 
        pprint(C)
        s = C.pop()
        if not any([g <= s for g in G]):
            if all(g != s for g in G):
                G.append(s)
            for g in G:
                vec = s.compute_s_vector(g)
                if vec and not any([g <= vec for g in G]):
                    C.append(vec)
                
            
                

    return [g for g in G if g[-1] >= 0]

