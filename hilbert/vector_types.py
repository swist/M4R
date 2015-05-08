from sympy import Matrix, pprint, floor
import itertools


class LiftableVector(Matrix):
    """docstring for LiftableVector"""
    linear_factors = []

    def __init__(self, *args, **kwargs):
        super(LiftableVector, self).__init__(*args, **kwargs)

    def __add__(self, other):
        z = super(LiftableVector, self).__add__(other)
        z.linear_factors = [x + y for x, y in itertools.izip(
            self.linear_factors,
            other.linear_factors
        )]
        return z

    def __sub__(self, other):
        z = super(LiftableVector, self).__sub__(other)
        z.linear_factors = [x - y for x, y in itertools.izip(
            self.linear_factors,
            other.linear_factors
        )]
        return z

    def __mul__(self, other):
        # if multiplying by a matrix, do whatever the superclass does
        # should not be needed
        if getattr(other, 'is_Matrix', False):
            return super(LiftableVector, self).__mul__(other)

        # otherwise multiply the linear factors as well
        z = super(LiftableVector, self).__mul__(other)
        z.linear_factors = [other * x for x in self.linear_factors]
        return z

    def __getitem__(self, key):
        from_super_class = super(LiftableVector, self).__getitem__(key)
        if(isinstance(key, slice)):
            if key.start or key.stop:
                return self.__class__(from_super_class)
        return from_super_class

    def lift(self, M):
        if M.rows == M.cols:
            i = 0
            h = self.col_join(Matrix([0]))
            while True:
                h[-1] = i
                echelon_form = M.row_join(h).rref()[0]
                coeffs = echelon_form[:,-1].values()
                if all(val.is_integer for val in coeffs):
                    return h
                i = i + 1
        else:
            return M * M[:M.cols,:].solve(self[:M.cols,:])
            

class BasisElement(LiftableVector):
    """docstring for BasisElement"""
    def __init__(self, *args, **kwargs):
        super(BasisElement, self).__init__(*args, **kwargs)

    def __le__(self, other):
        return all(
            [u_i <= v_i for u_i, v_i in itertools.izip(self[:-1], other[:-1])]
            ) and abs(self[-1]) <= abs(other[-1]) and self[-1]*other[-1] >= 0

    def __ge__(self, other):
        return other.__le__(self)

    def compute_s_vector(self, other):
        if self[-1]*other[-1] < 0:
            return self + other
        else:
            return None

    def compute_alpha(self, vector_g):
        return min(
            [floor(s_i/g_i) for s_i, g_i
                in zip(self[:], vector_g[:]) if g_i != 0]
        )

    def norm(self):
        return super(BasisElement, self).norm(1)

class ExtremeRay(LiftableVector):
    """docstring for ExtremeRay"""
    def __init__(self, *args, **kwargs):
        super(ExtremeRay, self).__init__(*args, **kwargs)

    def __le__(self, other):
        return self.support <= other.support

    @property
    def support(self):
        support = set()
        for index, el in enumerate(self):
            if el != 0:
                support.add(index)
        return support

    def norm(self):
        return len(self.support)

    def compute_s_vector(self, other):
        if self[-1]*other[-1] < 0:
            z = self - (self[-1]/other[-1])*other
            assert z[-1] == 0
            return z
        else:
            return

    def compute_alpha(self, vector_g):
        alpha = min([
            s_i/g_i for s_i, g_i in
            zip(self[0:-1], vector_g[0:-1]) if g_i != 0])
        assert alpha != 0
        return alpha

    def normal_form(self, set_G):
        s_now = self
        while True:
            if any([g == s_now for g in set_G]):
                return s_now - s_now

            normalizable = [g for g in set_G if (g <= s_now)]

            if not normalizable:
                break
            else:
                vector_g = normalizable[0]
                alpha = s_now.compute_alpha(vector_g)
                s_now = s_now - alpha * vector_g

        return s_now
