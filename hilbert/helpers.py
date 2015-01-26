import numpy as np
import math
import itertools
from sympy import *
from .basis_element import BasisElement
init_printing(use_unicode=True)

def compare_first_j_elements(u, v, j):
    #this can be parallelised
    return all([u_i <= v_i for u_i, v_i in itertools.izip(u[:j], v[:j])])

def compare_norms_of_jth_element(u, v, j):
    return abs(u[j]) <= abs(v[j])

def check_product_of_jth_elements_positive(u, v, j):
    return u[j]*v[j] >= 0

def check_vectors_square_leq(u, v, j):
    print "check_vectors_square_leq:"
    print "u:"
    pprint(u)
    print "v"
    pprint(v)
    print "j:"
    pprint(j)

    # these can be parallelised perhaps (run them asynchronously)
    return all([
        compare_first_j_elements(u, v, j),
        compare_norms_of_jth_element(u,v, j),
        check_product_of_jth_elements_positive(u, v, j)
        ])

def compute_alpha(vector_s, vector_g, j):
    print "compute_alpha:"
    print "s:"
    pprint(vector_s)
    print "v"
    pprint(vector_g)
    print "j:"
    pprint(j)

    return min([floor(s_i/g_i) for s_i, g_i in zip(vector_s[:j+2], vector_g[:j+2]) if g_i != 0])

def normal_form(vector_s, set_G, j):
    print "computing normal form from:"
    print "s:"
    pprint(vector_s)
    print "G:"
    pprint(set_G)

    while True:
        normalizable = [g for g in set_G if check_vectors_square_leq(g, vector_s, j) is True]
        print "normalizable"
        pprint(normalizable)
        if not normalizable:
            break
        else:
            for vector_g in normalizable:
                print "reducing s:"
                alpha = compute_alpha(vector_s, vector_g, j)
                print "alpha"
                print alpha
                vector_s = vector_s - alpha* vector_g

    return vector_s

def compute_s_vectors(vec_1, vec_2):
    if vec_1[-1]*vec_2[-1] < 0:
        return vec_1 + vec_2
    else:
        return None


def construct_c(G):
    C = []
    print "constructing C:\n"
    # use itertools product instead of zip for better performance
    for vector_f, vector_g in itertools.product(G, G):
        s_vector = compute_s_vectors(vector_f, vector_g)
        # only add if we don't have a vector that is equal
        if s_vector and not s_vector.is_zero and all((s_vector != c for c in C)):
            C.append(s_vector)
    print "C is:\n"
    pprint(C)
    return C

def extend_c(f,C,G):
    print "f not zero, extending C and G:\n"
    for g in G:
        vec = compute_s_vectors(f,g)
        print "computed s-vector:"
        pprint(vec)
        if vec and all((vec != c for c in C)):
            C.append(vec)

    if all(g != f for g in G):
        G.append(f)


def poittier(F, j):
    print "reduction algorithm:"
    G = F
    C = construct_c(G)
    print "reduction loop:"
    while len(C):
        print "C is:"
        pprint(C)
        print "G is:"
        pprint(G)
        s = C.pop()
        print "s is:"
        pprint(s)

        f = normal_form(s, G, j)

        print "f is:"
        pprint(f)
        if f and f.is_zero is False:
            extend_c(f,C,G)
        print "end of iteration \n"
    # return H^+
    return [g for g in G if g[-1] >= 0]

def ref(matrix):
    rows, cols = matrix.shape
    M = matrix.as_mutable()
    BC = ones(cols)
    m = min(cols, rows)
    i, stopper = 0,0
    while i < m  and stopper < cols:
        j = next((k for k in xrange(i, rows) if M[k,i] != 0), None)
        if j is not None:
            for k in xrange(j+1, rows):
                if M[k, i] != 0:
                    a = M[j,i]
                    b = M[k,i]
                    L = gcdex(a,b)
                    row_j = L[0]*M[j,:] + L[1]*M[k,:]
                    row_k = -b/L[2] * M[j,:] + a/L[2] * M[k,:]
                    M[j,:] = row_j
                    M[k,:] = row_k
            if i != j:
                M.row_swap(i, j)
            if M[i,i] < 0:
                M.row_op(i, lambda x, _: -1*x)
        else:
            M.col_swap(i, cols-1)
            BC.col_swap(i, cols-1)
            i = i - 1;
        i = i + 1 
        stopper = stopper + 1
    return matrix._new(M)