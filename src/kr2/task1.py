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
