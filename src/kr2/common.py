from sympy import *


class SplineData:
    slots = ("f", "xs")

    def __init__(self, f, xs):
        self.f = f
        self.xs = xs


def equispaced_nodes(a, b, n):
    return [a + i * (b - a) / (n - 1) for i in range(n)]


def chebyshev_nodes(a, b, n):
    return [
        (a + b) / 2 + (b - a) * cos((2 * k - 1) * pi / (2 * n)) / 2
        for k in range(1, n + 1)
    ]
