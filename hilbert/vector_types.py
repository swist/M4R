from sympy import ImmutableMatrix, pprint, floor
import itertools


class BaseVector(ImmutableMatrix):
    """docstring for BaseVector"""

    def __init__(self, *args, **kwargs):
        super(BaseVector, self).__init__(*args, **kwargs)

    def __getitem__(self, key):
        from_super_class = super(BaseVector, self).__getitem__(key)
        if(isinstance(key, slice)):
            if key.start or key.stop:
                return self.__class__(from_super_class)
        return from_super_class

    @property
    def origin(self):
        return self.as_mutable().is_zero

class BasisElement(BaseVector):
    """docstring for BasisElement"""
    def __init__(self, *args, **kwargs):
        super(BasisElement, self).__init__(*args, **kwargs)

    def __le__(self, other):
        return all(
            [u_i <= v_i for u_i, v_i in itertools.izip(self[:-1], other[:-1])]
            ) and abs(self[-1]) <= abs(other[-1]) and self[-1]*other[-1] >= 0

    def __ge__(self, other):
        return other.__le__(self)

    def s_vector(self, other):
        if self[-1]*other[-1] < 0:
            return self + other
        else:
            return BasisElement([0])

    def compute_alpha(self, vector_g):
        return min(
            [floor(s_i/g_i) for s_i, g_i
                in zip(self[:], vector_g[:]) if g_i != 0]
        )

    def lift(self, M):
        if M.rows == M.cols:
            i = 0
            h = self.as_mutable().col_join(ImmutableMatrix([0]))
            while True:
                h[-1] = i
                echelon_form = M.row_join(h).rref()[0]
                coeffs = echelon_form[:,-1].values()
                if all(val.is_integer for val in coeffs):
                    return BasisElement(h)
                i = i + 1
        else:
            return M * M[:M.cols,:M.cols].solve(self[:M.cols,:M.cols])            

    def norm(self):
        return super(BasisElement, self).norm(1)


    def normal_form(self, G):
        if any([g <= self for g in G]):
            return BasisElement([0])
        return self


class ExtremeRay(BaseVector):
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

    def lift(self, M):
        if M.rows == M.cols:
            comb = M[:self.rows, :self.rows].solve(self)
            return M[:, :self.rows] * comb

        else:
            return M * M[:M.cols,:M.cols].solve(self[:M.cols,:M.cols])  
    def norm(self):
        return len(self.support)

    def s_vector(self, other):
        if self[-1]*other[-1] < 0:
            z = self - (self[-1]/other[-1])*other
            assert z[-1] == 0
            return z
        else:
            return ExtremeRay([0])

    def compute_alpha(self, vector_g):
        alpha = min([
            s_i/g_i for s_i, g_i in
            zip(self[0:-1], vector_g[0:-1]) if g_i != 0])
        assert alpha != 0
        return alpha

    def normal_form(self, G):
        s = self
        print 'normal form extreme ray'
        pprint(G)
        pprint(s)
        while any(v <= s for v in G):
            g = next(v for v in G if v <= s)
            print(g, s)
            alpha = s.compute_alpha(g)
            s = s - alpha * g
            pprint(s)
        return s
