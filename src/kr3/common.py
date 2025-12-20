from sympy import *
from sympy.abc import x
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
    # return (
    #     sum(
    #         map(
    #             lambda k: f(a + (2 * k - 2) * h)
    #             + 4 * f(a + (2 * k - 1) * h)
    #             + f(a + 2 * k * h),
    #             range(1, n // 2 + 1),
    #         )
    #     )
    #     * h
    #     / 3
    # )
    res = 0

    for k in range(1, n + 1):
        res += f(a + (k - 1) * h) + 4 * f(a + (k - 0.5) * h) + f(a + k * h)

    return res * h / 6


def interpolation_quadrature(data: IntegralData):
    a, b, f, n = data.a, data.b, data.f, data.n
    h = (b - a) / n
    lagrange_data = LagrangeData([a + k * h for k in range(1, n + 1)], f)
    p = []

    for k in range(n):
        basis = lagrange_basis_poly(lagrange_data, k)
        p.append(integrate(basis, (x, a, b)))

    res = 0

    for k in range(1, n + 1):
        res += f(a + k * h) * p[k - 1]

    return res


# def gauss_quadrature(data: IntegralData):
#     a, b, f, n = data.a, data.b, data.f, data.n
#     h = (b - a) / n
