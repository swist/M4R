from .. import helpers
from ..vector_types import BasisElement


def test_construct_C():
    G = [
        BasisElement([[1, 13]]),
        BasisElement([[0, 18]]),
        BasisElement([[0, -18]])
    ]

    C = [
        BasisElement([[1, -5]])
    ]

    result = helpers.construct_c(G)

    assert result == C


def test_extend_c():
    f = BasisElement([[1, -5]])
    C = []
    G = [
        BasisElement([[1, 13]]),
        BasisElement([[0, 18]]),
        BasisElement([[0, -18]])
    ]

    Gnew = [
        BasisElement([[1, 13]]),
        BasisElement([[0, 18]]),
        BasisElement([[0, -18]]),
        f
    ]

    Cnew = [
        BasisElement([[2, 8]]),
        BasisElement([[1, 13]])
    ]

    helpers.extend_c(f, C, G)
    assert C == Cnew
    assert G == Gnew


def test_normal_form():
    G = [
        BasisElement([[1, 13]]),
        BasisElement([[0, 18]]),
        BasisElement([[0, -18]])
    ]

    s = BasisElement([[1, -5]])

    f = BasisElement([[1, -5]])

    result = helpers.normal_form(s, G, 1)
    assert f == result


def test_poittier():
    F = [
        BasisElement([[1, 13]]),
        BasisElement([[0, 18]]),
        BasisElement([[0, -18]])
    ]

    result = [
        BasisElement([[1, 13]]),
        BasisElement([[0, 18]]),
        BasisElement([[2, 8]]),
        BasisElement([[3, 3]]),
        BasisElement([[7, 1]]),
        BasisElement([[18, 0]])
    ]

    H2 = helpers.poittier(F, 1)
    assert H2 == result
