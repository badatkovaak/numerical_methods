from sympy import *
from sympy.abc import x
from src.kr2.task1 import create_first_order_spline
from src.kr2.task2 import create_third_order_spline
from src.kr2.common import *


if __name__ == "__main__":
    a, b, n = 0, 1, 11
    f = Lambda(x, ((x - 1) / 3) ** 8)
    data = SplineData(f, equispaced_nodes(a, b, n))
    g = create_first_order_spline(data)
    h = create_third_order_spline(data)
    # plot(f, g, h, (x, a, b))
