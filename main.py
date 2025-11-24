from sympy import *
from sympy.abc import x
from src.kr2.task1 import first_order_spline
from src.kr2.task2 import (
    third_order_spline_by_definition,
    third_order_spline_by_moments,
)
from src.kr2.common import *


def task1():
    a, b, n = 0, 1, 5
    f = Lambda(x, ((x - 1) / 2) ** 9)

    data_equispaced = SplineData(f, equispaced_nodes(a, b, n))
    data_chebyshev = SplineData(f, chebyshev_nodes(a, b, n))

    spline_equispaced = first_order_spline(data_equispaced)
    spline_chebyshev = first_order_spline(data_chebyshev)
    print(spline_equispaced, "\n")
    print(spline_chebyshev, "\n")
    print()

    control_points = [0.29, 0.42, 0.76]
    diffs_equispaced = compare_at_control_points(
        f, Lambda(x, spline_equispaced), control_points
    )
    diffs_chebyshev = compare_at_control_points(
        f, Lambda(x, spline_chebyshev), control_points
    )
    print(diffs_equispaced)
    print(diffs_chebyshev)
    print()


def task2():
    a, b, n = 0, 1, 5
    f = Lambda(x, ((x - 1) / 2) ** 9)
    g = Lambda(x, x**2)
    data = SplineData(g, equispaced_nodes(a, b, n))
    spline = third_order_spline_by_moments(data)
    # plot(g, Lambda(x, spline), (x, a, b))


if __name__ == "__main__":
    task2()
