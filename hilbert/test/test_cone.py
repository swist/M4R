from ..cone import Cone
from sympy import Matrix, pprint
from ..vector_types import BasisElement


def test_return_rays():
    C = Cone([[1, 2], [1, -3]])
    rays = C.rays
    assert rays == [C.row(0), C.row(1)]


def test_dual():
    C = Cone([[1, 2], [1, -3]])
    dual = C.dual

    pprint(dual)
    assert dual == Cone([[2, -1], [3, 1]])

    assert C.dual.dual == C

    C = Cone([
        [1, 1, 0],
        [0, 3, 0],
        [0, 0, 1]
    ])

    assert C.dual == Cone([
        [1, 0, 0],
        [-1, 1, 0],
        [0, 0, 1]
    ])

    C = Cone([[1, 2], [2, 1]])
    assert C.dual == Cone([[2, -1], [-1, 2]])


def test_points_contain():
    C = Cone([[1, 2], [2, 1]])
    v = Matrix([1, 0])
    assert C.contains(v) == False

    v = Matrix([2, 0])
    pprint(v.normalized())
    assert C.contains(v) == False
    v = Matrix([4, 2])
    assert C.contains(v) == True

    v = Matrix([3, 0])
    assert C.contains(v) == False


    C = Cone([[1, 2], [0, 3]])
    v = Matrix([1, 0])
    assert C.contains(v) == False
    v = Matrix([1, 1])
    assert C.contains(v) == False
    v = Matrix([1, 2])
    assert C.contains(v) == True
    assert C.contains(v.T) == True



def test_dual_big():
    C = Cone([
        [0, 0, 0, 0, 0, 3, -4, -1, 2],
        [0, 0, 0, 0, 1, -1, 1, 0, -1],
        [0, 0, 0, 1, 2, 0, 0, -1, -2],
        [0, 0, 1, 0, 1, 0, 0, -1, -1],
        [0, 1, 2, 0, 0, 0, 0, -1, -2],
        [1, 0, 2, 0, 0, 0, 0, -2, -1],
        [-2, 0, -2, 0, 0, 0, 0, 3, 0],
        [-2, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, -2, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, -1, 0]
    ])
    assert len(C.dual.rays) == 10
    assert C.dual == Cone([
        [0, -2,  0,  0, -1, -2, -2,  0, -1],
        [0,  1,  0,  0,  0,  0,  0,  0,  0],
        [0,  0,  0, -2,  1,  4,  3,  0,  0],
        [0,  0,  0,  1,  0,  0,  0,  0,  0],
        [0,  0,  0,  0,  0, -4, -3,  0,  0],
        [0,  0,  0,  0,  0, -1, -1,  0,  0],
        [-1,  0, -2, -2, -1,  0,  0, -2, -1],
        [-2,  0, -1,  0, -1, -2, -1, -2,  0],
        [-1, -2,  0,  0, -1, -2, -2,  0, -1],
        [0, -2, -1, -2, -1,  0, -1,  0, -2]
        ])


def test_hilbert_basis():
    C = Cone([[1, 2], [1, -3]])

    basis = C.hilbert_basis()
    pprint(basis)

    assert basis == [
        BasisElement([[1, 1]]),
        BasisElement([[1, 2]]),
        BasisElement([[1, 0]]),
        BasisElement([[1, -1]]),
        BasisElement([[1, -2]]),
        BasisElement([[1, -3]])
    ]


def test_hilbert_basis2():
    C = Cone([[1, 2], [2, 1]])

    basis = C.hilbert_basis()
    pprint(basis)

    assert basis == [
        BasisElement([[1, 1]]),
        BasisElement([[1, 2]]),
        BasisElement([[2, 1]]),
    ]


def test_hilbert_basis3():
    C = Cone([[1, 0], [0, 1]])

    basis = C.hilbert_basis()
    pprint(basis)

    assert basis == [
        BasisElement([[1, 0]]),
        BasisElement([[0, 1]])
    ]


def test_hilbert_basis4():
    C = Cone([[2, 1], [1, -1]])

    basis = C.hilbert_basis()
    pprint(basis)

    assert basis == [
        BasisElement([[1, 0]]),
        BasisElement([[2, 1]]),
        BasisElement([[1, -1]])
    ]


def test_hilbert_basis_3d2():
    C = Cone([
        [1, 1, 0],
        [0, 3, 0],
        [0, 0, 1]
    ])
    basis = C.hilbert_basis()
    pprint(basis)
    result = [
        BasisElement([[1, 1, 0]]),
        BasisElement([[0, 1, 0]]),
        BasisElement([[0, 0, 1]])
    ]
    assert len(basis) == len(result)
    assert basis == result


def test_hilbert_basis_3d():
    C = Cone([
        [1, 1, 0],
        [0, 3, 0],
        [0, 0, 1],
        [2, 5, 11]
    ])
    basis = C.hilbert_basis()
    pprint(basis)
    result = [
        BasisElement([[1, 1, 0]]),
        BasisElement([[0, 1, 0]]),
        BasisElement([[0, 0, 1]])
    ]
    assert len(basis) == len(result)
    assert basis == result


def test_hilbert_basis_3d3():
    C = Cone([
        [1, 1, 0],
        [0, 3, 0],
        [0, 0, 1],
        [2, -5, 11]
    ])

    assert C.dual == Cone([
            [1,  0, 0],
            [0,  0, 1],
            [-11, 11, 7],
            [5,  2, 0]
        ])

    basis = C.hilbert_basis()
    pprint(basis)
    result = [
        BasisElement([[1,  1,  0]]),
        BasisElement([[0,  1,  0]]),
        BasisElement([[0,  0,  1]]),
        BasisElement([[2, -5, 11]]),
        BasisElement([[1,  0,  2]]),
        BasisElement([[1, -2,  5]]),
        BasisElement([[2, -3,  8]]),
        BasisElement([[1, -1,  4]])
    ]
    assert len(basis) == len(result)
    assert basis == result


# def test_hilbert_basis_big():
#     C = Cone([
#         [0, 0, 0, 0, 0, 3, -4, -1, 2],
#         [0, 0, 0, 0, 1, -1, 1, 0, -1],
#         [0, 0, 0, 1, 2, 0, 0, -1, -2],
#         [0, 0, 1, 0, 1, 0, 0, -1, -1],
#         [0, 1, 2, 0, 0, 0, 0, -1, -2],
#         [1, 0, 2, 0, 0, 0, 0, -2, -1],
#         [-2, 0, -2, 0, 0, 0, 0, 3, 0],
#         [-2, 0, 0, 0, 0, 0, 0, 1, 0],
#         [0, 0, -2, 0, 0, 0, 0, 1, 0],
#         [0, 0, 0, 0, 0, 0, 0, -1, 0]
#     ])
#     pprint(C.rows)
#     basis = C.hilbert_basis()
#     pprint(basis)

#     result = [
#         BasisElement([1, 0, 2, 2, 1, 0, 0, 2, 1]),
#         BasisElement([2, 0, 1, 0, 1, 2, 1, 2, 0]),
#         BasisElement([1, 2, 0, 0, 1, 2, 2, 0, 1]),
#         BasisElement([1, 1, 1, 1, 1, 1, 1, 1, 1]),
#         BasisElement([0, 2, 1, 2, 1, 0, 1, 0, 2])
#     ]

#     assert basis == result
