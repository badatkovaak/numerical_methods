from sympy import *
from math import fabs

x = symbols("x")


class LagrangeData:
    __slots__ = ("points", "f")

    def __init__(self, xs, f):
        self.points = xs
        self.f = f


def lagrange_basis_poly(data: LagrangeData, k: int):
    res = 1

    for i, x_i in enumerate(data.points):
        if i == k:
            continue

        res *= (x - x_i) / (data.points[k] - x_i)

    return res


def lagrange_poly(data: LagrangeData):
    result = 0

    for i, x_i in enumerate(data.points):
        result += data.f(x_i) * lagrange_basis_poly(data, i)

    return simplify(result).collect(x)


def equidistant_nodes():
    return [k / 4 for k in range(5)]


def chebyshev_nodes():
    return [1 / 2 + cos((2 * k - 1) * pi / 10) / 2 for k in range(1, 6)]


def compute_values_at_control_points(data: LagrangeData, lagrange_poly):
    control_points = [0.29, 0.42, 0.76]

    for i, x_i in enumerate(control_points):
        value_of_poly = lagrange_poly.as_poly().eval(x_i).as_expr().round(15)
        value_of_func = data.f(x_i)
        diff = fabs(value_of_poly - value_of_func)

        print(f"for {i + 1} control point:")
        print(f"\tvalue of lagrange polynomial is {value_of_poly}")
        print(f"\tvalue of the function is {value_of_func}")
        print(f"\tabsolute difference is {diff}")


def compute_max_omega(data: LagrangeData):
    omega = Poly(1, x)

    for x_i in data.points:
        omega *= x - x_i

    omega_prime = Poly(omega.diff().as_expr())
    roots = list(
        filter(
            lambda y: y >= 0 and y <= 1,
            map(lambda x: x.evalf(14), omega_prime.real_roots()),
        )
    )
    roots.append(0)
    roots.append(1)
    max_omega = 0

    for root in roots:
        if max_omega < fabs(omega.eval(root)):
            max_omega = omega.eval(root)

    return max_omega


if __name__ == "__main__":
    f = Lambda(x, ((x - 1) / 3) ** 8)

    data_equi = LagrangeData(equidistant_nodes(), f)
    lagrange_poly_equi = lagrange_poly(data_equi)
    print(lagrange_poly_equi)
    compute_values_at_control_points(data_equi, lagrange_poly_equi)
    print(f"maximum of omega is {compute_max_omega(data_equi)}")

    print()

    data_chebyshev = LagrangeData(chebyshev_nodes(), f)
    lagrange_poly_chebyshev = lagrange_poly(data_chebyshev)
    print(lagrange_poly_chebyshev)
    compute_values_at_control_points(data_chebyshev, lagrange_poly_chebyshev)

    print()

    print(f"maximum of omega is {compute_max_omega(data_chebyshev)}")
