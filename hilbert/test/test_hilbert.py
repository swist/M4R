from sympy import Matrix, FiniteSet, init_printing, pprint
from ..bases_computation import construct_hilbert_basis
from ..helpers import ref
from ..basis_element import BasisElement

init_printing()


# def test_hilbert_basis_construction1():
#     result = [
#         BasisElement([[1,0,2,2,1,0,0,2,1]]),
#         BasisElement([[2,0,1,0,1,2,1,2,0]]),
#         BasisElement([[1,2,0,0,1,2,2,0,1]]),
#         BasisElement([[1,1,1,1,1,1,1,1,1]]),
#         BasisElement([[0,2,1,2,1,0,1,0,2]])
#     ]

    # M = Matrix([
    #     [0,0,0,0,0,3,-4,-1,2],
    #     [0,0,0,0,1,-1,1,0,-1],
    #     [0,0,0,1,2,0,0,-1,-2],
    #     [0,0,1,0,1,0,0,-1,-1],
    #     [0,1,2,0,0,0,0,-1,-2],
    #     [1,0,2,0,0,0,0,-2,-1],
    #     [-2,0,-2,0,0,0,0,3,0],
    #     [-2,0,0,0,0,0,0,1,0],
    #     [0,0,-2,0,0,0,0,1,0],
    #     [0,0,0,0,0,0,0,-1,0]
    # ])

#     basis = construct_hilbert_basis(M)
#     pprint(basis)
#     pprint(result)

#     assert len(basis) == len(result)

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

    assert len(basis) == len(result)

def test_hilbert_basis_construction2():
    M = Matrix([[1,2], [1, -3]])
    result = [
        BasisElement([[1,-3]]),
        BasisElement([[1,-2]]),
        BasisElement([[1,-1]]),
        BasisElement([[1,0]]),
        BasisElement([[1,1]]),
        BasisElement([[1,2]])
    ]

    basis = construct_hilbert_basis(M)
    pprint(basis)
    pprint(result)

    assert len(basis) == len(result)


def test_hilbert_basis_construction3():
    M = Matrix([
        [1,1,0],
        [0,3,0],
        [0,0,1]
    ]);
    basis = construct_hilbert_basis(M)
    result = [        
        BasisElement([1,1,0]),
        BasisElement([0,1,0]),        
        BasisElement([0,0,1])
    ]
    assert len(basis) == len(result)