from sympy import *
from typing import List
from math import fabs

class LagrangeData:
    __slots__ = ('points', 'f', 'x')
    
    def __init__(self, xs: List[float], f):
        self.points = xs
        self.f = f
        self.x = symbols('x')

def create_lagrange_basis_poly(data: LagrangeData, k: int):
    res = 1
    x = data.x

    for (i, x_i) in enumerate(data.points):
        if i == k:
            continue
        
        res *= (x - x_i)/(data.points[k] - x_i)

    return res.expand()

def create_lagrange_poly(data: LagrangeData):
    result = 0

    for (i, x_i) in enumerate(data.points):
        result += data.f(x_i) * create_lagrange_basis_poly(data, i)

    return result.expand()

def equidistant_nodes():
    return [0, 1/4, 1/2, 3/4, 1]

def chebyshev_nodes():
    return [1/2 + cos((2*k - 1)*pi/10)/2 for k in range(1, 6)]

def compute_at_control_points(data: LagrangeData, lagrange_poly):
    control_points = [0.29, 0.42, 0.76]

    for (i, x_i) in enumerate(control_points):
        value_of_poly = lagrange_poly.eval(x_i)
        value_of_func = data.f(x_i)
        diff = fabs(value_of_poly - value_of_func)

        print(f"for {i + 1} control point:")
        print(f"\tvalue of lagrange polynomial is {value_of_poly}")
        print(f"\tvalue of the function is {value_of_func}")
        print(f"\tabsolute difference is {diff}")

def finite_difference(data: LagrangeData, m: int, k: int):
    result = 0

    for j in range(m):
        result += (-1)**j * binomial(m, j)*data.f(data.points[k + m - j])

    return result

def omega(data: LagrangeData, k: int):
    result = 1

    for i in range(k):
        pass


if __name__ == "__main__":
    # x = symbol('x')
    # data = LagrangeData([1,2,3], Lambda(x:= symbols('x'), ((x-1)/3)**8))
    data = LagrangeData([1,2,3], Lambda(x:= symbols('x'), x**2))
    print(create_lagrange_poly(data))
    print(chebyshev_nodes())

    # for node in chebyshev_nodes():
    #     print(node >=0, node <= 1)
    pass
