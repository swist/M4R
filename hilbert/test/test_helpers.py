from sympy import Matrix, pprint
from .. import helpers
from ..vector_types import BasisElement

def test_ref_2d():
    M = Matrix([[1,2], [2,1]])
    result, pivots = helpers.ref(M)

    assert result == Matrix([[1,2], [0, 3]])
    assert pivots == Matrix.eye(2)
    

def test_ref_3d():
    M_3= Matrix([
        [1, 2, 3, 1], 
        [2, 3, 1, 1], 
        [3, 1, 2, 1]
    ])
    result_3, pivots_3 = helpers.ref(M_3)

    assert result_3 == Matrix([
        [1, 2, 3, 1], 
        [0, 1, 5, 1], 
        [0, 0, 18, 3]
    ])
    assert pivots_3 == Matrix.eye(4)




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


def test_poittier():
    F = [
        BasisElement([[1,0]]),
        BasisElement([[0,5]]),
        BasisElement([[0,-5]])
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

    pprint(H2)
    assert H2 == result



