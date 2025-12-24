from sympy import *
from sympy.abc import x, y
from sympy.printing.latex import latex
from functools import reduce
from src.kr1.task1 import *


class IntegralData:
    __slots__ = ("f", "a", "b", "n")

    def __init__(self, f, a: float, b: float, n: int):
        self.f = f
        self.a = a
        self.b = b
        self.n = n


def left_rectangles(data: IntegralData):
    h = (data.b - data.a) / data.n
    return sum(map(lambda k: data.f(data.a + k * h), range(data.n))) * h


def right_rectangles(data: IntegralData):
    h = (data.b - data.a) / data.n
    return sum(map(lambda k: data.f(data.a + (k + 1) * h), range(data.n))) * h


def middle_rectangles(data: IntegralData):
    h = (data.b - data.a) / data.n
    return sum(map(lambda k: data.f(data.a + (k + 0.5) * h), range(data.n))) * h


def simpsons_formula(data: IntegralData):
    a, f, n = data.a, data.f, data.n
    h = (data.b - data.a) / n
    res = 0

    for k in range(1, n + 1):
        res += f(a + (k - 1) * h) + 4 * f(a + (k - 0.5) * h) + f(a + k * h)

    return res * h / 6


def interpolation_quadrature(data: IntegralData, points=None, ro=None):
    a, b, f, n = data.a, data.b, data.f, data.n

    if points is None:
        h = (b - a) / n
        points = [a + k * h for k in range(n)]

    ro = ro if ro is not None else 1
    lagrange_data = LagrangeData(points, f)
    p = []

    for k in range(n):
        basis = lagrange_basis_poly(lagrange_data, k)
        p.append(integrate(basis * ro, (x, a, b)))

    return sum(map(lambda k: f(points[k]) * p[k], range(n))).evalf(14)


def gauss_quadrature(data: IntegralData, omega, ro):
    a, b = data.a, data.b

    if a != -1 or b != 1:
        e = data.f.expr.subs(x, x * (b - a) / 2 + (a + b) / 2)
        data.f = Lambda(x, e)
        data.a, data.b = -1, 1

    omega = Poly(omega.as_expr())
    roots = list(map(lambda x: x.evalf(14), omega.real_roots()))

    return interpolation_quadrature(data, roots, ro).evalf(14)


expr1 = sqrt(1.2 * x + 0.7) / (1.4 * x + sqrt(1.3 * x * x + 0.5))
expr2 = sin(0.7 * x + 0.4) / (2.2 + cos(0.3 * x * x + 0.7))

a1, b1, f1 = (
    1.2,
    2.8,
    Lambda(x, expr1),
)
a2, b2, f2 = (
    0.5,
    1.3,
    Lambda(x, expr2),
)


def task1():
    n = 10
    data1 = IntegralData(f1, a1, b1, n)
    data2 = IntegralData(f2, a2, b2, n)

    print()
    print(left_rectangles(data1))
    print(right_rectangles(data1))
    print()
    print(left_rectangles(data2))
    print(right_rectangles(data2))
    print()


def task2():
    n = 10
    data1 = IntegralData(f1, a1, b1, n)
    data2 = IntegralData(f2, a2, b2, n)

    print()
    print(middle_rectangles(data1))
    print()
    print(middle_rectangles(data2))
    print()


def task3():
    n = 10
    data1 = IntegralData(f1, a1, b1, n)

    print()
    print(simpsons_formula(data1))
    print()


def task4():
    n = 6
    data2 = IntegralData(f2, a2, b2, n)
    omega = legendre(n, x)
    ro = 1

    print()
    print(gauss_quadrature(data2, omega, ro))
    print()


def task5():
    n = 10
    data2 = IntegralData(f2, a2, b2, n)

    print()
    print(interpolation_quadrature(data2))
    print()


def task6():
    n = 6
    data1 = IntegralData(f1, a1, b1, n)
    ro = 1 / sqrt(1 - x * x)
    omega = chebyshevt(n, x)

    print()
    print(gauss_quadrature(data1, omega, ro))
    print()
