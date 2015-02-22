from sympy import Matrix, pprint
from hilbert.bases_computation import construct_hilbert_basis
from hilbert.vector_types import ExtremeRay


class Cone(Matrix):
    """docstring for Cone"""
    def __init__(self, *args, **kwargs):
        super(Cone, self).__init__(*args, **kwargs)

    def rays(self):
        return [self.row(i) for i in range(0, self.rows)]

    def contains(self, vector):
        v = vector.normalized()
        for col in range(0, self.cols):
            vals = self.col(col).values()
            vals.append(0)
            if v[col] < min(vals):
                return False
            if v[col] > max(vals):
                return False
        return True

    def dual(self):
        rays = construct_hilbert_basis(self.T, ExtremeRay)
        return [self.row_join(ray.T).nullspace()[0] for ray in rays]
