from sympy import ImmutableMatrix, pprint, eye, gcd, Matrix
from hilbert.bases_computation import construct_generating_set
from hilbert.hnf import hnf_row, hnf_col
from hilbert.vector_types import ExtremeRay, BasisElement


def preimage(matrix, vectors, VectorClass):
    nullspaces = [(-1 * v).row_join(matrix).nullspace() for v in vectors]
    
    nullspaces = [to_postive_halfspace(ns[0]) for ns in nullspaces]
    vectors = [VectorClass(ns[1:matrix.cols+1]) for ns in nullspaces]
    pprint(vectors)
    return set(map(convert_to_ZZ, vectors))

def to_postive_halfspace(vector):
    # make sure 1 in the first position
    if vector[0] < 0:
        vector = vector * abs(vector[0])/vector[0]
    return vector

def convert_to_ZZ(vector):
    pprint(vector)
    for i in range(len(vector)):
        if(vector[i].is_Rational):
            vector = vector[i].q * vector
    pprint(vector)
    return vector


class Cone(ImmutableMatrix):
    """docstring for Cone"""

    _dual = None
    _rays = None

    def __init__(self, *args, **kwargs):
        super(Cone, self).__init__(*args, **kwargs)

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
            self._dual = Cone(col_matrix)
            self._dual._dual = self
        return self._dual
        

    def _compute_dual(self):  
        pprint(self)   
        A, BC = self.rref()
        pprint(A)
        dual_rays = construct_generating_set(Matrix(A), ExtremeRay)
        dual_rays = preimage(self.T, dual_rays, ExtremeRay)

        return dual_rays

    def hilbert_basis(self):
        dual = self._compute_dual()
        return self.hb_from_dual(dual)


    def hb_from_dual(self, dual):
        A, BC = hnf_row(dual.T)
        pprint(A)
        pprint(BC)
        basis = construct_generating_set(A)
        print('basis')
        pprint(basis)
        basis = preimage(dual, basis, BasisElement)
        pprint(basis)
        return basis

   