import pytest
from ..cone import Cone
from sympy import Matrix, pprint
from ..vector_types import BasisElement, ExtremeRay


def test_return_rays():
    C = Cone([[1, 1], [2, -3]])
    rays = C.rays
    assert rays == [C.col(0), C.col(1)]

def test_generators_in_hnf():
    C = Cone([[1, 1], [2, -3]])
    hnf = C.to_hnf()
    assert hnf == Matrix([[1,1],[0,5]])

    C = Cone([
        [1, 1, 0],
        [0, 3, 0],
        [0, 0, 1]
    ])

    four_gens = Matrix([
        [1, 1, 0],
        [0, 3, 0],
        [0, 0, 1],
        [2, -5, 11],
    ]).T

    C = Cone(four_gens)
    assert C.to_hnf() == Matrix([
        [1, 3, 0, -5],
        [0, 3, 0, -7],
        [0, 0, 1, 11]
    ])

def test_hb_from_dual2D():
    # symmetric 2D cone
    dual_in = Matrix([
        [2, -1],
        [-1, 2]
    ])

    C = Cone(dual_in)

    res_exp = [
        BasisElement([1,2]),
        BasisElement([2,1]),
        BasisElement([1,1])
    ]
    res_bas = C.hb_from_dual(dual_in)
    assert len(res_exp) == len(res_bas)
    assert set(map(tuple, res_exp)) == set(map(tuple,res_bas))

    dual_in = Matrix([
        [2,-1],
        [3, 1]
    ])

    # another 2D cone with elts in negative orthant
    C = Cone(dual_in)

    res_exp = [
        BasisElement([1,-3]),
        BasisElement([1,-2]),
        BasisElement([1,-1]),
        BasisElement([1,0]),
        BasisElement([1,1]),
        BasisElement([1,2])
    ]
    res_bas = C.hb_from_dual(dual_in)
    
    assert len(res_exp) == len(res_bas)
    assert set(map(tuple, res_exp)) == set(map(tuple,res_bas))


    dual_in = Matrix([
        [1, 0], 
        [0, 1]
    ])
    C = Cone(dual_in)

    res_bas = C.hb_from_dual(dual_in)   

    res_exp = [
        BasisElement([1, 0]),
        BasisElement([0, 1])
    ]

    assert len(res_exp) == len(res_bas)
    assert set(map(tuple, res_exp)) == set(map(tuple,res_bas))

def test_hilbert_basis_big():
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
    dual_in = Matrix([
        [ 0, -2,  0,  0, -1, -2, -2,  0, -1],
        [ 0,  1,  0,  0,  0,  0,  0,  0,  0],
        [ 0,  0,  0, -2,  1,  4,  3,  0,  0],
        [ 0,  0,  0,  1,  0,  0,  0,  0,  0],
        [ 0,  0,  0,  0,  0, -4, -3,  0,  0],
        [ 0,  0,  0,  0,  0, -1, -1,  0,  0],
        [-1,  0, -2, -2, -1,  0,  0, -2, -1],
        [-2,  0, -1,  0, -1, -2, -1, -2,  0],
        [-1, -2,  0,  0, -1, -2, -2,  0, -1],
        [ 0, -2, -1, -2, -1,  0, -1,  0, -2]
    ])

    res_exp = [
        BasisElement([ 0, 0,  0, 0, 0,  3, -4, -1,  2]),
        BasisElement([ 0, 0,  0, 0, 1, -1,  1,  0, -1]),
        BasisElement([ 0, 0,  0, 1, 2,  0,  0, -1, -2]),
        BasisElement([ 0, 0,  1, 0, 1,  0,  0, -1, -1]),
        BasisElement([ 0, 1,  2, 0, 0,  0,  0, -1, -2]),
        BasisElement([ 1, 0,  2, 0, 0,  0,  0, -2, -1]),
        BasisElement([-2, 0, -2, 0, 0,  0,  0,  3,  0]),
        BasisElement([-2, 0,  0, 0, 0,  0,  0,  1,  0]),
        BasisElement([ 0, 0, -2, 0, 0,  0,  0,  1,  0]),
        BasisElement([ 0, 0,  0, 0, 0,  0,  0, -1,  0]),
        BasisElement([-1, 0, -1, 0, 0,  0,  0,  1,  0]),
        BasisElement([-1, 0, -2, 0, 0,  0,  0,  2,  0]),
        BasisElement([-2, 0, -1, 0, 0,  0,  0,  2,  0]),
        BasisElement([ 0, 0, -1, 0, 0,  0,  0,  0,  0]),
        BasisElement([-1, 0,  0, 0, 0,  0,  0,  0,  0])
    ]
    res_bas = C.hb_from_dual(dual_in)
    assert len(res_exp) == len(res_bas)
    assert set(map(tuple, res_exp)) == set(map(tuple,res_bas))

def test_hilbert_basis_3d3():
    C = Cone([
        [1, 1, 0],
        [0, 3, 0],
        [0, 0, 1],
        [2, -5, 11]
    ])

    dual_in = Matrix([
        [1,  0, 0],
        [0,  0, 1],
        [-11, 11, 7],
        [5,  2, 0]
    ])

    res_bas = C.hb_from_dual(dual_in)
    
    res_exp = [
        BasisElement([[1,  1,  0]]),
        BasisElement([[0,  1,  0]]),
        BasisElement([[0,  0,  1]]),
        BasisElement([[2, -5, 11]]),
        BasisElement([[1,  0,  2]]),
        BasisElement([[1, -2,  5]]),
        BasisElement([[2, -3,  8]]),
        BasisElement([[1, -1,  4]])
    ]
    assert len(res_exp) == len(res_bas)
    assert set(map(tuple, res_exp)) == set(map(tuple,res_bas))
    assert False


def test_rank_gt_1():
    dual_in = Matrix([
        [13, 0, 1], 
        [23, 1, 0], 
        [0, 0, 1], 
        [-11, 11, 7], 
        [5, 2, 0]
    ])
    C = Cone(dual_in)
    res_bas = C.hb_from_dual(dual_in)

    res_exp = [
        BasisElement([1, 1, 0]),
        BasisElement([0, 1, 0]),
        BasisElement([0, 0, 1]),
        BasisElement([2, -5, 11]),
        BasisElement([-1, 23, 13]),
        BasisElement([2, -3, 8]),
        BasisElement([1, 0, 2]),
        BasisElement([1, -2, 5]),
        BasisElement([1, -1, 4])
    ]
    assert len(res_exp) == len(res_bas)
    assert set(map(tuple, res_exp)) == set(map(tuple,res_bas))


@pytest.mark.broken
def test_dual():
    C = Cone([[1,2],[2,1]])
    dual = C._compute_dual()    
    dual_rays = [
        ExtremeRay([2,-1]),
        ExtremeRay([-1,2]),
    ]
    assert set(map(tuple, dual)) == set(map(tuple, dual_rays))

    C = Cone([[1, 1], [2, -3]])
    dual = C._compute_dual()
    pprint(dual)
    dual_rays = [
        ExtremeRay([2,-1]),
        ExtremeRay([3,1]),
    ]
    assert set(map(tuple, dual)) == set(map(tuple, dual_rays))

    M = Matrix([
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [2, -5, 11]
    ]).T

    C = Cone(M)

    dual = C._compute_dual()
    pprint(dual)
    dual_rays = [
        ExtremeRay([1, 0, 0]),
        ExtremeRay([0, 0, 1]),
        ExtremeRay([-11, 11, 7]),
        ExtremeRay([5, 2, 0])
    ]
    assert set(map(tuple, dual)) == set(map(tuple, dual_rays))

    C = Cone([
        [1, 0, 0],
        [1, 3, 0],
        [0, 0, 1]
    ])
    dual = C._compute_dual()
    dual_rays = [
        ExtremeRay([1, 0, 0]),
        ExtremeRay([-1, 1, 0]),
        ExtremeRay([0, 0, 1])
    ]
    print('dual')
    pprint(dual)
    assert set(map(tuple, dual)) == set(map(tuple, dual_rays))

@pytest.mark.broken
def test_dual_big():
    C = Cone(Matrix([
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
    ]).T)
    
    dual = C._compute_dual()
    
    dual_rays = [
        ExtremeRay([0, -2,  0,  0, -1, -2, -2,  0, -1]),
        ExtremeRay([0,  1,  0,  0,  0,  0,  0,  0,  0]),
        ExtremeRay([0,  0,  0, -2,  1,  4,  3,  0,  0]),
        ExtremeRay([0,  0,  0,  1,  0,  0,  0,  0,  0]),
        ExtremeRay([0,  0,  0,  0,  0, -4, -3,  0,  0]),
        ExtremeRay([0,  0,  0,  0,  0, -1, -1,  0,  0]),
        ExtremeRay([-1,  0, -2, -2, -1,  0,  0, -2, -1]),
        ExtremeRay([-2,  0, -1,  0, -1, -2, -1, -2,  0]),
        ExtremeRay([-1, -2,  0,  0, -1, -2, -2,  0, -1]),
        ExtremeRay([0, -2, -1, -2, -1,  0, -1,  0, -2])
    ]
    assert set(map(tuple, dual)) == set(map(tuple, dual_rays))




