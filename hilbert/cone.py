from sympy import ImmutableMatrix, pprint, eye, gcd
from hilbert.bases_computation import construct_generating_set
from hilbert.hnf import hnf_row, hnf_col
from hilbert.fm import fourierMotzkin
from hilbert.vector_types import ExtremeRay, BasisElement


def preimage(matrix, vectors, VectorClass):
    nullspaces = [(-1 * v).row_join(matrix).nullspace() for v in vectors]
    
    nullspaces = [convert_to_ZZ(ns[0]) for ns in nullspaces]
    vectors = [VectorClass(ns[1:matrix.cols+1]) for ns in nullspaces]
    pprint(vectors)
    return vectors


def convert_to_ZZ(vector):
    for i in range(len(vector)):
        if(vector[i].is_Rational):
            vector = vector[i].q * vector
    # make sure 1 in the first position
    if vector[0] < 0:
        vector = vector * abs(vector[0])/vector[0]
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
            pprint([gcd(*self.col(i)) for i in range(self.cols)])
            self._rays = [ExtremeRay(self.col(i)) for i in range(self.cols)]
        return self._rays

    def is_dual_ray(self, dual_ray):
        rays = self.rays
        return any([dual_ray.dot(ray) == 0 for ray in rays])

    def to_hnf(self):
        return hnf_row(self)[0]

    def _compute_dual(self):
        if self._dual is None:
            # pprint(self)
            # dual_rays = fourierMotzkin(self)

            A, BC = hnf_row(self)
            pprint(A)
            pprint(BC)
            if self.rows is 2:
                dr = A.inv() * BC.inv()
                self._dual = map(convert_to_ZZ, [dr.row(i) for i in range(dr.rows)])
            else:
                dual_rays = construct_generating_set(A, ExtremeRay)
                dual_rays = preimage(self.T, dual_rays, ExtremeRay)
                pprint(dual_rays)
                self._dual = dual_rays


        return self._dual


    def hilbert_basis(self):
        dual = self._compute_dual()
        return self.hb_from_dual(d)


    def hb_from_dual(self, dual):
        A, BC = hnf_row(dual.T)
        pprint(A)
        pprint(BC)
        basis = construct_generating_set(A, BasisElement)
        print('basis')
        pprint(basis)
        basis = preimage(dual, basis, BasisElement)
        return basis

   