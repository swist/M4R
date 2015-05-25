from ..vector_types import BasisElement, LiftableVector
from sympy import Matrix, Rational

def test_basis_leq():
    g = BasisElement([[1, 13]])
    s = BasisElement([[1, 13]])

    assert g[:2] <= s[:2]
    s = BasisElement([[1, 14]])

    assert g[:2] <= s[:2]

    s = BasisElement([[1, -14]])
    assert not (g[:2] <= s[:2])

    g = BasisElement([[0, 1, 13]])
    s = BasisElement([[1, 0, 3]])
    assert not (g[:3] <= s[:3])


def test_lift():
    s = LiftableVector([1])
    m = Matrix([[1, 2], [0, 3]])
    result = LiftableVector([1, 0])

    lift = s.lift(m)
    assert lift == result

    M = Matrix([
        [1,0,0],
        [5,2,0],
        [-11,11,7]
    ])
    H = [
        LiftableVector([1, 1]),
        LiftableVector([0, 2]),
        LiftableVector([2, 0]),
    ]

    res = [
        LiftableVector([1, 1, 2]),
        LiftableVector([0, 2, 4]),
        LiftableVector([2, 0, 0]),
    ]

    lifted = [h.lift(M) for h in H]
    assert res == lifted





def test_lift_multiple_choice():
    M = Matrix([
        [1, 5, -11, 0],
        [0, 2, 11, 0],
        [0, 0, 7, 1]
    ]).T

    H = [
        LiftableVector([1,1,2]),
        LiftableVector([0,2,4]),
        LiftableVector([2,0,0]),
        LiftableVector([0,0,7]),
        LiftableVector([0,4,1]),
        LiftableVector([1,7,0]),
        LiftableVector([0,14,0])
    ]

    exp = [
        LiftableVector([1,1,2,5]),
        LiftableVector([0,2,4,-1]),
        LiftableVector([2,0,0,11]),
        LiftableVector([0,0,7,1]),
        LiftableVector([0,4,1,-3]),
        LiftableVector([1,7,0,0]),
        LiftableVector([0,14,0,-11])
    ]


    lifted = [h.lift(M) for h in H]
    assert set(map(tuple,exp)) == set(map(tuple, exp))



def test_lift_casting():    
    s = LiftableVector([1])
    m = Matrix([[1, 2], [0, 3]])
    result = LiftableVector([1, 0])

    lift = s.lift(m)
    assert lift == result
    assert type(lift) == LiftableVector

    
    s = BasisElement([1])
    m = Matrix([[1, 2], [0, 3]])
    result = BasisElement([1, 0])
    lift = s.lift(m)
    assert lift == result
    assert type(lift) == BasisElement
