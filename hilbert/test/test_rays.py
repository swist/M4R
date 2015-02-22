from sympy import pprint
from ..cone import Cone


def test_dual_rays():
    C = Cone([
        [0, 1],
        [3, 2]
    ])

    rays = C.dual()

    pprint(rays)
    assert False
