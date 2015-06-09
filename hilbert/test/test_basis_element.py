from ..vector_types import BasisElement
from sympy import Matrix, Rational

def test_basis_leq():
    g = BasisElement([[1, 13]])
    s = BasisElement([[1, 13]])

    assert g <= s
    s = BasisElement([[1, 14]])

    assert g[:2] <= s[:2]

    s = BasisElement([[1, -14]])
    assert not (g[:2] <= s[:2])

    g = BasisElement([[0, 1, 13]])
    s = BasisElement([[1, 0, 3]])
    assert not g <= s

def test_addition():
    g = BasisElement([[1, 13]])
    s = BasisElement([[1, 23]])
    assert g+s == BasisElement([[2, 36]])


def test_s_vector():
    g = BasisElement([[1, -1]])
    s = BasisElement([[1, 2]])

    assert g.s_vector(s) == BasisElement([[2,1]])

def test_lift():
    s = BasisElement([1])
    m = Matrix([[1, 2], [0, 3]])
    result = BasisElement([1, 0])

    lift = s.lift(m)
    assert lift == result

    M = Matrix([
        [1,0,0],
        [5,2,0],
        [-11,11,7]
    ])
    H = {
        BasisElement([1, 1]),
        BasisElement([0, 2]),
        BasisElement([2, 0]),
    }

    res = {
        BasisElement([1, 1, 2]),
        BasisElement([0, 2, 4]),
        BasisElement([2, 0, 0]),
    }

    lifted = {h.lift(M) for h in H}
    assert res == lifted





def test_lift_multiple_choice():
    M = Matrix([
        [1, 5, -11, 0],
        [0, 2, 11, 0],
        [0, 0, 7, 1]
    ]).T

    H = {
        BasisElement([1,1,2]),
        BasisElement([0,2,4]),
        BasisElement([2,0,0]),
        BasisElement([0,0,7]),
        BasisElement([0,4,1]),
        BasisElement([1,7,0]),
        BasisElement([0,14,0])
    }

    exp = {
        BasisElement([1,1,2,5]),
        BasisElement([0,2,4,-1]),
        BasisElement([2,0,0,11]),
        BasisElement([0,0,7,1]),
        BasisElement([0,4,1,-3]),
        BasisElement([1,7,0,0]),
        BasisElement([0,14,0,-11])
    }


    lifted = {h.lift(M) for h in H}
    assert exp == lifted

def test_is_origin():
    f = BasisElement([0])
    assert f.origin is True


def test_lift_casting():    
    s = BasisElement([1])
    m = Matrix([[1, 2], [0, 3]])
    result = BasisElement([1, 0])

    lift = s.lift(m)
    assert lift == result
    assert type(lift) == BasisElement

    
    s = BasisElement([1])
    m = Matrix([[1, 2], [0, 3]])
    result = BasisElement([1, 0])
    lift = s.lift(m)
    assert lift == result
    assert type(lift) == BasisElement
