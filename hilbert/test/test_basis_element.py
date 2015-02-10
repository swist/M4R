from ..vector_types import BasisElement, ExtremeRay
from sympy import Matrix, Rational

refed_matrix = Matrix([
    [1, 2, 3],
    [0, 1, -1143],
    [0, 0, 6921]
])


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


def test_ray_leq():
    g = ExtremeRay([[1, 0]])
    s = ExtremeRay([[1, 13]])

    assert g[:2] <= s[:2]

    g = ExtremeRay([[1, 0]])
    assert g[:2] <= s[:2]

    g = ExtremeRay([[0, 1, 13]])
    s = ExtremeRay([[1, 0, 3]])
    assert not (g[:3] <= s[:3])

    g = ExtremeRay([[1, 0, 0]])
    assert g[:3] <= s[:3]


def test_ray_support():
    v = ExtremeRay([[1, 2, 3]])
    assert v.support() == {0, 1, 2}

    v = ExtremeRay([[0, 1, -1]])
    assert v.support() == {1, 2}


def test_ray_s_vector():
    v = ExtremeRay([[0, -5]])
    w = ExtremeRay([[1, 2]])

    s = v.compute_s_vector(w)
    # always 0 in the last place
    assert s[-1] == 0


def test_ray_alpha():
    g = ExtremeRay([[1, 5, 132]])
    s = ExtremeRay([[1, 3, 13]])

    assert s.compute_alpha(g) == Rational(3, 5)

    s = ExtremeRay([[1, 0]])
    g = ExtremeRay([[Rational(5, 2), 1]])

    assert s.compute_alpha(g) == Rational(2, 5)


def test_constructor():
    h = BasisElement([refed_matrix[0, 0]])
    assert h.cols == 1
    assert h.rows == 1
    # test if it constructed itself with only one linear factor
    assert len(h.linear_factors) == 1
    # test if it only has one linear factor for the first element
    assert h.linear_factors[0] is h.linear_factors[-1]
    # check if it takes 1 as the default linear factor
    assert h.linear_factors[-1] == 1

    h_2 = BasisElement([refed_matrix[1, :2]])
    assert h_2.cols == 2
    assert h_2.rows == 1
    assert len(h_2.linear_factors) == 2
    assert h_2.linear_factors[0] is not h_2.linear_factors[-1]
    assert h_2.linear_factors[-1] == 1
    assert h_2.linear_factors[0] == 0

    assert type(h[:2]) is BasisElement
    assert h[:2] == h


def test_addition():
    h = BasisElement([refed_matrix[0, :2]])
    h_2 = BasisElement([refed_matrix[1, :2]])

    z = h + h_2

    assert type(z) is BasisElement
    assert z.linear_factors[0] == h.linear_factors[0] + h_2.linear_factors[0]
    assert z.linear_factors[1] == h.linear_factors[1] + h_2.linear_factors[1]
    assert len(z.linear_factors) == len(h.linear_factors) == len(h_2.linear_factors)


def test_subtraction():
    h = BasisElement([refed_matrix[0, :2]])
    h_2 = BasisElement([refed_matrix[1, :2]])

    z = h - h_2

    assert type(z) is BasisElement
    assert z.linear_factors[0] == h.linear_factors[0] - h_2.linear_factors[0]
    assert z.linear_factors[1] == h.linear_factors[1] - h_2.linear_factors[1]
    assert len(z.linear_factors) == len(h.linear_factors) == len(h_2.linear_factors)


def test_scalar_multiplication():
    h = BasisElement([refed_matrix[0, :2]])
    z = 3 * h

    assert type(z) is BasisElement
    assert z == 3*h
    assert z.linear_factors[0] == 3 * h.linear_factors[0]
    assert z.linear_factors[1] == 3 * h.linear_factors[1]
    assert len(z.linear_factors) == len(h.linear_factors)


def test_lift_multiple_choice():
    h = BasisElement([refed_matrix[0, 0]])

    h = h.lift_multiple_choice(refed_matrix[:2, :2])
    assert h.cols == 2
    assert h[0, -1] == 0
    assert h.linear_factors[-1] == -2
    h = h.lift_multiple_choice(refed_matrix[:3, :3])
    assert h.cols == 3
    assert h[0, -1] == 2289
    assert h.linear_factors[-1] == 0

    h = BasisElement([[1]])
    h = h.lift_multiple_choice(Matrix([[0], [1]]))

    assert h.cols == 2
    assert h[0, -1] == 0
    assert h.linear_factors[-1] == 0
