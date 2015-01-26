from sympy import Matrix, FiniteSet, init_printing, pprint
from ..bases_computation import construct_hilbert_basis
from ..helpers import ref
from ..basis_element import BasisElement

init_printing()


def test_hilbert_basis_construction1():
    result = [
        BasisElement([[1,0,2,2,1,0,0,2,1]]),
        BasisElement([[2,0,1,0,1,2,1,2,0]]),
        BasisElement([[1,2,0,0,1,2,2,0,1]]),
        BasisElement([[1,1,1,1,1,1,1,1,1]]),
        BasisElement([[0,2,1,2,1,0,1,0,2]])
    ]

    M = Matrix([
        [0,0,0,0,0,3,-4,-1,2],
        [0,0,0,0,1,-1,1,0,-1],
        [0,0,0,1,2,0,0,-1,-2],
        [0,0,1,0,1,0,0,-1,-1],
        [0,1,2,0,0,0,0,-1,-2],
        [1,0,2,0,0,0,0,-2,-1],
        [-2,0,-2,0,0,0,0,3,0],
        [-2,0,0,0,0,0,0,1,0],
        [0,0,-2,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,-1,0]
    ])

    basis = construct_hilbert_basis(M)
    pprint(basis)
    pprint(result)

    assert basis == result

def test_hilbert_basis_construction():
    result = [
        BasisElement([[1,1]]),
        BasisElement([[0,3]]),
        BasisElement([[3,0]])
    ]

    M = Matrix([[1,-2], [-2,1]])

    basis = construct_hilbert_basis(M)
    pprint(basis)
    pprint(result)

    assert basis == result

# def test_lifting():
#   K_j_plus_one = Matrix([[1,-2], [0, 3]])
#   H = FiniteSet(Matrix([[1]]))

#   result = lift(H,K_j_plus_one)
#   assert result == FiniteSet()
