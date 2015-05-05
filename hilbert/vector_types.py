from sympy import Matrix, pprint, floor
import itertools


class LiftableVector(Matrix):
    """docstring for LiftableVector"""
    linear_factors = []

    def __init__(self, *args, **kwargs):
        super(LiftableVector, self).__init__(*args, **kwargs)
        self.linear_factors = [0] * (self.cols - 1)
        self.linear_factors.append(1)

    def __getitem__(self, key):
        from_super_class = super(LiftableVector, self).__getitem__(key)
        if(isinstance(key, slice)):
            if key.start or key.stop:
                return self.__class__(from_super_class)
        return from_super_class

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

    def lift_single_choice(self, M):
        i = 0
        h = self.col_join(Matrix([0]))
        while True:
            h[-1] = i
            echelon_form = M.row_join(h).rref()[0]
            coeffs = echelon_form[:,-1].values()
            if all(val.is_integer for val in coeffs):
                return h
            i = i + 1


    def lift_multiple_choice(self, row_to_add):
        self.linear_factors.append(0)
        linear_factors_vector = Matrix([self.linear_factors])

        last_el = (linear_factors_vector * row_to_add)[-1]
        while last_el < 0:
            last_el = last_el + row_to_add[-1]
            self.linear_factors[-1] = self.linear_factors[-1] + 1

        if last_el >= row_to_add[-1]:

            self.linear_factors[-1] = -1 * (last_el/row_to_add[-1])
            last_el = last_el % row_to_add[-1]
        # crazy mutation stuff ahppening here
        linear_factors = self.linear_factors
        self = self.row_join(Matrix([last_el]))
        self.linear_factors = linear_factors

        return self


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
        return sum(self.values())


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

    def compute_s_vector(self, other):
        if self[-1]*other[-1] < 0:
            z = self - (self[-1]/other[-1])*other
            return z
        else:
            return

    def lift(self, A):
        return self.lift_single_choice(A)

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
