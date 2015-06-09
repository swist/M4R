from sympy import ImmutableMatrix, pprint, eye, gcd, Matrix
from hilbert.bases_computation import (
    construct_generating_set, preimage, to_postive_halfspace, convert_to_ZZ) 
from hilbert.hnf import hnf_row, hnf_col
from hilbert.vector_types import ExtremeRay, BasisElement

class Cone(ImmutableMatrix):
    """docstring for Cone"""

    _dual = None
    _rays = None

    def __init__(self, *args, **kwargs):
        super(Cone, self).__init__()

    @property
    def rays(self):
        if self._rays is None:
            self._rays = {ExtremeRay([self.col(i)]) for i in range(self.cols)}
        return self._rays

    @property 
    def dual(self):
        if self._dual is None:
            rays = self._compute_dual()
            col_matrix = Matrix.hstack(*[ray.as_mutable() for ray in rays])
            self._dual = Cone(col_matrix.T)
        return self._dual
        

    def _compute_dual(self):  
        pprint(self)   
        A, BC = self.rref()
        pprint(A)
        dual_rays = construct_generating_set(Matrix(A), ExtremeRay)
        dual_rays = preimage(self.T, dual_rays, ExtremeRay)

        return dual_rays

    def hilbert_basis(self):
        dual = self.dual
        return self._hb_from_dual(dual)


    def _hb_from_dual(self, dual):
        A, BC = hnf_row(dual.T)
        pprint(A)
        pprint(BC)
        basis = construct_generating_set(A)
        print('basis')
        pprint(basis)
        basis = preimage(dual, basis, BasisElement)
        pprint(basis)
        return basis

   