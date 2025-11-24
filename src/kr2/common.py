from sympy import *
from typing import Callable
from math import fabs


class SplineData:
    slots = ("f", "xs")

    def __init__(self, f: Callable[[float], float], xs: list[float]):
        self.f = f
        self.xs = xs


def equispaced_nodes(a: float, b: float, n: int):
    return [a + i * (b - a) / (n - 1) for i in range(n)]


def chebyshev_nodes(a: float, b: float, n: int):
    return [
        (a + b) / 2 + (b - a) * cos((2 * k - 1) * pi / (2 * n)) / 2
        for k in range(1, n + 1)
    ]


def compare_at_control_points(
    f: Callable[[float], float],
    g: Callable[[float], float],
    control_points: list[float],
) -> list[float]:
    diffs = []

    for point in control_points:
        diffs.append(fabs(f(point) - g(point)))

    return diffs
