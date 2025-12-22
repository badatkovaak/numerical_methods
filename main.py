from sympy import *
from sympy.abc import x
from src.kr3.common import *


if __name__ == "__main__":
    expr = x * x
    f, a, b, n = Lambda(x, expr), -1, 1, 3
    # omega = legendre(n, x)
    omega = chebyshevt(n, x)
    ro = 1 / sqrt(1 - x * x)
    print(omega)
    data = IntegralData(f, a, b, n)
    print(integrate(expr * ro, (x, a, b)).evalf(14))
    print()
    # print(left_rectangles(data))
    # print(right_rectangles(data))
    # print(middle_rectangles(data))
    # print(simpsons_formula(data))
    # print(interpolation_quadrature(data))
    print(gauss_quadrature(data, omega, ro))
