from sympy import Matrix
from .. import helpers
from ..basis_element import BasisElement

def test_ref_2d():
    M = Matrix([[1,-2], [-2,1]])
    result = Matrix([[1,-2], [0, 3]])

    assert helpers.ref(M) == result
    assert not M is result

def test_ref_3d():
    M_3= Matrix([
        [1, 2, 3, 1], 
        [2, 3, 1, 1], 
        [3, 1, 2, 1]
    ])
    result_3= Matrix([
        [1, 2, 3, 1], 
        [0, 1, 5, 1], 
        [0, 0, 18, 3]
    ])

    assert helpers.ref(M_3) == result_3
    assert not M_3 is result_3


def test_leq_helpers():
    V_1 = Matrix([0,1,3]).T
    V_2 = Matrix([0,2,0]).T

    assert helpers.compare_first_j_elements(V_1, V_2, 2) == True
    assert helpers.check_product_of_jth_elements_positive(V_1, V_2, 2) == True

    V_2 = Matrix([1,0,3]).T
    assert helpers.compare_first_j_elements(V_1, V_2, 2) == False


def test_construct_C():
    G = [
        BasisElement([[1,13]]),
        BasisElement([[0,18]]),
        BasisElement([[0,-18]])
    ]

    C = [
        BasisElement([[1,-5]])
    ]

    result = helpers.construct_c(G)

    assert result == C

def test_extend_c():
    f = BasisElement([[1,-5]])
    C = []
    G = [
        BasisElement([[1,13]]),
        BasisElement([[0,18]]),
        BasisElement([[0,-18]])
    ]

    Gnew = [
        BasisElement([[1,13]]),
        BasisElement([[0,18]]),
        BasisElement([[0,-18]]),
        f
    ]

    Cnew = [
        BasisElement([[2,8]]),
        BasisElement([[1,13]])
    ]

    helpers.extend_c(f, C, G)
    assert C == Cnew
    assert G == Gnew

def test_normal_form():
    G = [
        BasisElement([[1,13]]),
        BasisElement([[0,18]]),
        BasisElement([[0,-18]])
    ]

    s = BasisElement([[1,-5]])

    f = BasisElement([[1,-5]])

    result = helpers.normal_form(s,G,1)
    assert f == result


def test_poittier():
    F = [
        BasisElement([[1,13]]),
        BasisElement([[0,18]]),
        BasisElement([[0,-18]])
    ]

    result = [
        BasisElement([[1,13]]),
        BasisElement([[0,18]]),
        BasisElement([[2,8]]),
        BasisElement([[3,3]]),
        BasisElement([[7,1]]),
        BasisElement([[18,0]])
    ]

    H2 = helpers.poittier(F,1)
    assert H2 == result



