from sympy import ImmutableMatrix, pprint
from hilbert.bases_computation import construct_generating_set
from hilbert.hnf import hnf_row
from hilbert.vector_types import ExtremeRay, BasisElement


def preimage(matrix, vectors, VectorClass):
    if matrix.rows != vectors[0].rows:
        vectors = [v.T for v in vectors]

    nullspaces = [matrix.row_join(-1 * v).nullspace() for v in vectors]
    nullspaces = [ns for ns in nullspaces if ns]
    print('nsses')
    print([len(ns) for ns in nullspaces])
    vectors = [VectorClass([ns[0][:matrix.cols]]) for ns in nullspaces]
    print(vectors)
    return vectors


def convert_to_ZZ(vector):
    for i in range(len(vector)):
        if(vector[i].is_Rational):
            vector = vector[i].q * vector
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
            self._rays = [ExtremeRay(self.row(i)) for i in range(0, self.rows)]

        return self._rays

    @property
    def dual(self):
        if self._dual is None:
            self._compute_dual()
        return self._dual

    def contains(self, vector):

        dual_rays = self.dual.rays
        pprint(vector)
        pprint(dual_rays)
        pprint([vector.dot(ray) >= 0 for ray in dual_rays])
        return all([vector.dot(ray) >= 0 for ray in dual_rays])

    def is_dual_ray(self, dual_ray):
        rays = self.rays
        return any([dual_ray.dot(ray) == 0 for ray in rays])

    def _compute_dual(self):
        A, BC = hnf_row(self.T)
        print(BC)

        rays = construct_generating_set(A, ExtremeRay)
        pprint(rays)

        rays = preimage(self, rays, ExtremeRay)

        rays = map(convert_to_ZZ, rays)
        pprint(rays)
        rays = [ray for ray in rays if self.is_dual_ray(ray)]
        print('rays')
        pprint(rays)
        self._dual = Cone(rays)
        self._dual._dual = self
        return self._dual

    def hilbert_basis(self):
        dual = self.dual
        A, BC = hnf_row(dual.T)

        h_n = construct_generating_set(A, BasisElement)
        print(h_n)
        print('basis el len')
        pprint(len(h_n[0]))
        nullspaces = [dual.row_join(-1 * v.T).nullspace() for v in h_n]
        print(nullspaces)
        rays = [ns[0][:self.cols] for ns in nullspaces if ns]
        rs = [BasisElement([ray]) for ray in rays]
        print('basis')
        pprint(rs)
        return rs
