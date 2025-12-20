from sympy import *
from sympy.abc import x
from src.kr2.task1 import first_order_spline
from src.kr2.task2 import (
    third_order_spline_by_definition,
    third_order_spline_by_moments,
)
from src.kr2.common import *
from src.kr3.common import *


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


def task2():
    a, b, n = 0, 1, 5
    f = Lambda(x, ((x - 2) / 2) ** 9)

    data_equispaced = SplineData(f, equispaced_nodes(a, b, n))
    data_chebyshev = SplineData(f, chebyshev_nodes(a, b, n))

    print(spline_equispaced, spline_chebyshev, sep="\n\n", end="\n\n")

    control_points = [0.29, 0.42, 0.76]
    diffs_equispaced = compare_at_control_points(
        f, Lambda(x, spline_equispaced), control_points
    )
    diffs_chebyshev = compare_at_control_points(
        f, Lambda(x, spline_chebyshev), control_points
    )
    print(diffs_equispaced, diffs_chebyshev, sep="\n", end="\n\n")


def compute_maxes():
    a, b, n = 0, 1, 5

    print(compute_max_omega(equispaced_nodes(a, b, n)))
    print(compute_max_omega(chebyshev_nodes(a, b, n)))


if __name__ == "__main__":
    # task1()
    # task2()
    # compute_maxes()
    f, a, b = Lambda(x, x * x * x * x * x), 1, 3
    data = IntegralData(f, a, b, 6)
    print(integrate(x * x * x * x * x, (x, a, b)))
    # print(left_rectangles(data))
    # print(right_rectangles(data))
    # print(middle_rectangles(data))
    # print(simpsons_formula(data))
    print(interpolation_quadrature(data))
