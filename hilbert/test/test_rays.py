from sympy import pprint, Rational, Matrix
from ..vector_types import ExtremeRay
from ..cone import Cone


def test_le():
    v1 = ExtremeRay([1, -1])
    v2 = ExtremeRay([1, 0])
    v3 = ExtremeRay([0, -3])

    assert v1 >= v2

    v1 = v1 - v2
    pprint(v1.support)
    pprint(v2.support)

    pprint(v1.support & v2.support)
    assert not v1 <= v2
    assert v1 <= v3


def test_ray_support():
    v = ExtremeRay([[1, 2, 3]])
    assert v.support == {0, 1, 2}

    v = ExtremeRay([[0, 1, -1]])
    assert v.support == {1, 2}


def test_ray_s_vector():
    v = ExtremeRay([[0, -5]])
    w = ExtremeRay([[1, 2]])

    s = v.s_vector(w)
    # always 0 in the last place
    assert s[-1] == 0


def test_ray_alpha():
    g = ExtremeRay([[1, 5, 132]])
    s = ExtremeRay([[1, 3, 13]])

    assert s.compute_alpha(g) == Rational(3, 5)

    s = ExtremeRay([[1, 0]])
    g = ExtremeRay([[Rational(13, 132), 1]])

    assert s.compute_alpha(g) == Rational(132, 13)


def test_normal_form():
    G = [
        ExtremeRay([1, -1]),
        ExtremeRay([0, 3]),
        ExtremeRay([0, -3])
    ]

    s = ExtremeRay([1, 0])

    res = s.normal_form(G)
    assert res == s


def test_lift():
    H = {
        ExtremeRay([1,1]),
        ExtremeRay([0,5]),
    }

    K = Matrix([
        [1,1,0,0],
        [0,5,-1,-7],
        [0,0,2,3]
        ]).T

    lifted = {h.lift(K[:3,:3]) for h in H}
    expected = {ExtremeRay([1,1,0]), ExtremeRay([0,5,-1])}
    assert expected == lifted
