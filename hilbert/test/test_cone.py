from ..cone import Cone
from sympy import Matrix, pprint


def test_return_rays():
    C = Cone([[1, 2], [1, -3]])
    rays = C.rays()
    assert rays == [C.row(0), C.row(1)]


# def test_points_contain():
#   C = Cone([[1,2],[1,-3]])
#   v = Matrix([1,0])
#   assert C.contains(v) == True

#   v = Matrix([2,0])
#   pprint(v.normalized())
#   assert C.contains(v) == True
#   v = Matrix([1,-4])
#   assert C.contains(v) == False

def test_dual():
    C = Cone([[1, 2], [1, -3]])
    dual = C.dual()

    pprint(dual)
    assert C.dual == Cone([[2, -1], [3, 1]])

    assert C.dual().dual() == C
