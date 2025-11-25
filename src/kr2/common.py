from sympy import *
from sympy.abc import x
from typing import Callable
from math import fabs
from functools import reduce


class SplineData:
    slots = ("f", "xs")

    def __init__(self, f: Callable[[float], float], xs: list[float]):
        self.f = f
        self.xs = xs


def equispaced_nodes(a: float, b: float, n: int):
    return [a + i * (b - a) / (n - 1) for i in range(n)]


def chebyshev_nodes(a: float, b: float, n: int):
    return sorted(
        [
            (a + b) / 2 + (b - a) * cos((2 * k - 1) * pi / (2 * n)) / 2
            for k in range(1, n + 1)
        ]
    )


def compare_at_control_points(
    f: Callable[[float], float],
    g: Callable[[float], float],
    control_points: list[float],
) -> list[float]:
    return list(map(lambda p: fabs(f(p) - g(p)), control_points))


def compute_max_omega(points):
    omega = reduce(lambda z, p: z * (x - p), points, sympify(1))
    omega_prime = Poly(omega.diff().expand())
    roots = map(lambda z: z.evalf(), omega_prime.real_roots())
    return max(map(lambda z: abs(omega.subs(x, z).evalf()), roots))
