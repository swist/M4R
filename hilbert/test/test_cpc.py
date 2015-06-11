from ..vector_types import BasisElement, BaseVector
from ..bases_computation import construct_generating_set, critical_pairs
from sympy import Matrix, pprint

def test_construct_C():
    G = {
        BasisElement([[1, 13]]),
        BasisElement([[0, 18]]),
        BasisElement([[0, -18]])
    }

    C = {
        BasisElement([[1, -5]])
    }
    result = reduce(set.union, (critical_pairs(f, G) for f in G))

    assert result == C

def test_cpc():
	matrix_in = Matrix([
		[1, 5, -11, 0],
		[0, 2, 11, 0],
		[0, 0, 7, 1]
		])

	res_exp = {
		BasisElement([1,1,2,5]),
		BasisElement([2,0,0,11]),
		BasisElement([0,0,7,1]),
		BasisElement([1, 7, 0, 0]),
		BasisElement([2, 4, 1, 8]),
		BasisElement([1, 5, 3, 2]),
		BasisElement([1, 3, 6, 4]),
		BasisElement([0, 2, 11, 0])
	}

	res_bas = construct_generating_set(matrix_in)
	pprint(res_bas)
	
	assert len(res_exp) == len(res_bas)
	assert res_exp == res_bas

def test_HB_2D_hnf():
	matrix_in = Matrix([
		[1,-2],
		[0, 3]
	])

	res_exp = {
		BasisElement([1,1]),
		BasisElement([0,3]),
		BasisElement([3,0])
	}

	res_bas = construct_generating_set(matrix_in)
	pprint(res_bas)
	
	assert len(res_exp) == len(res_bas)
	assert res_exp == res_bas

	matrix_in = Matrix([
		[1, -1],
		[0, 5]
	])

	res_exp = {
		BasisElement([1,4]),
		BasisElement([0,5]),
		BasisElement([2,3]),
		BasisElement([3,2]),
		BasisElement([4,1]),
		BasisElement([5,0])
	}

	res_bas = construct_generating_set(matrix_in)
	pprint(res_bas)
	
	assert len(res_exp) == len(res_bas)
	assert res_exp == res_bas




