from sympy import *
from sympy.abc import x
from src.kr2.common import SplineData
from typing import Any


def first_order_base_spline(data: SplineData, k: int) -> Any:
    xs = data.xs

    if k == 0:
        cond1 = (0, x < xs[0])
        cond2 = ((xs[1] - x) / (xs[1] - xs[0]), x <= xs[1])
        cond3 = (0, True)
        return Piecewise(cond1, cond2, cond3)

    if k == (len(xs) - 1):
        cond1 = (0, x < xs[k - 1])
        cond2 = ((x - xs[k - 1]) / (xs[k] - xs[k - 1]), x <= xs[k])
        cond3 = (0, True)
        return Piecewise(cond1, cond2, cond3)

    cond1 = (0, x < xs[k - 1])
    cond2 = ((x - xs[k - 1]) / (xs[k] - xs[k - 1]), x <= xs[k])
    cond3 = ((xs[k + 1] - x) / (xs[k + 1] - xs[k]), x < xs[k + 1])
    cond4 = (0, True)
    return Piecewise(cond1, cond2, cond3, cond4)


def first_order_spline(data: SplineData) -> Any:
    result = 0

    for k, x_k in enumerate(data.xs):
        b = create_first_order_base_spline(data, k)
        result += data.f(x_k) * b

    return piecewise_exclusive(piecewise_fold(result))
