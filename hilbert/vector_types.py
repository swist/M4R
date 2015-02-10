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

    def support(self):
        support = set()
        for index, el in enumerate(self):
            if el != 0:
                support.add(index)
        return support

    def lift_single_choice(self, row_to_add):
        print("self")
        pprint(self)
        row_to_add = row_to_add[:, -1]
        print("linear factors")
        print(self.linear_factors)
        linear_factors_vector = Matrix([self.linear_factors])
        print("row to add")
        pprint(row_to_add)
        print("linear factos vec")
        pprint(linear_factors_vector)

        last_el = (linear_factors_vector * row_to_add)[-1]
        linear_factors = self.linear_factors
        self = self.row_join(Matrix([last_el]))
        self.linear_factors = linear_factors
        return self

    def lift_multiple_choice(self, row_to_add):
        self.linear_factors.append(0)
        linear_factors_vector = Matrix([self.linear_factors])
        print("row to add")
        pprint(row_to_add)
        last_el = (linear_factors_vector * row_to_add)[-1]
        while last_el < 0:
            last_el = last_el + row_to_add[-1]
            self.linear_factors[-1] = self.linear_factors[-1] + 1
        pprint(last_el)
        if last_el >= row_to_add[-1]:
            pprint("last >= row")
            self.linear_factors[-1] = -1 * (last_el/row_to_add[-1])
            last_el = last_el % row_to_add[-1]

        print("last_el")
        pprint(last_el)
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


class ExtremeRay(LiftableVector):
    """docstring for ExtremeRay"""
    def __init__(self, *args, **kwargs):
        super(ExtremeRay, self).__init__(*args, **kwargs)

    def __le__(self, other):
        supp_1 = self.support()
        supp_2 = other.support()
        print("lhs support")
        print(supp_1)
        print("rhs support")
        print(supp_2)
        return supp_1.issubset(supp_2)

    def compute_s_vector(self, other):
        if self[-1]*other[-1] < 0:
            z = self - (self[-1]/other[-1])*other
            return z
        else:
            return

    def compute_alpha(self, vector_g):
        pprint(self)
        pprint(vector_g)
        self_copy = self[:-1]
        g_copy = vector_g[:-1]
        pprint(self_copy)
        pprint(g_copy)
        return min([
            s_i/g_i for s_i, g_i in
            zip(self_copy[:], g_copy[:]) if g_i != 0])
