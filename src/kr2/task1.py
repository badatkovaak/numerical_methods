from sympy import *
from sympy.abc import x
from src.kr2.common import SplineData
from typing import Any


def first_order_spline(data: SplineData) -> Any:
    xs = data.xs
    f = data.f
    conds = [(0, x < xs[0])]

    for k in range(len(xs) - 1):
        cond = f(xs[k]) * (x - xs[k + 1]) / (xs[k] - xs[k + 1]) + f(xs[k + 1]) * (
            x - xs[k]
        ) / (xs[k + 1] - xs[k])
        conds.append((cond, x <= xs[k + 1]))

    conds.append((0, True))

    return Piecewise(*conds).evalf()


def task1():
    a, b, n = 0, 1, 5
    f = Lambda(x, ((x - 1) / 2) ** 9)

    data_equispaced = SplineData(f, equispaced_nodes(a, b, n))
    data_chebyshev = SplineData(f, chebyshev_nodes(a, b, n))

    spline_equispaced = first_order_spline(data_equispaced)
    spline_chebyshev = first_order_spline(data_chebyshev)
    print(spline_equispaced, spline_chebyshev, sep="\n\n", end="\n\n")

    control_points = [0.29, 0.42, 0.76]
    diffs_equispaced = compare_at_control_points(
        f, Lambda(x, spline_equispaced), control_points
    )
    diffs_chebyshev = compare_at_control_points(
        f, Lambda(x, spline_chebyshev), control_points
    )
    print(diffs_equispaced, diffs_chebyshev, sep="\n", end="\n\n")
