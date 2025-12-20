from sympy import *

x = symbols("x")


class LagrangeData:
    __slots__ = ("points", "f")

    def __init__(self, xs, f):
        self.points = xs
        self.f = f


def create_lagrange_basis_poly(data: LagrangeData, k: int):
    res = 1

    for i, x_i in enumerate(data.points):
        if i == k:
            continue

        res *= (x - x_i) / (data.points[k] - x_i)

    return res


class HermitData:
    __slots__ = ("points", "f", "f_prime")

    def __init__(self, points, f, f_prime=None):
        self.points = points
        self.f = f

        if f_prime is None:
            self.f_prime = Lambda(x, Derivative(f.expr, x, evaluate=True))
        else:
            self.f_prime = f_prime


def omega(data: HermitData):
    result = 1

    for x_i in data.points:
        result *= x - x_i

    return result


def compute_c(data: HermitData, k: int):
    return omega(data).diff().diff().as_poly().eval(data.points[k]) / omega(
        data
    ).diff().as_poly().eval(data.points[k])


def create_hermit_poly(data: HermitData):
    result = 0

    for i, x_i in enumerate(data.points):
        l = create_lagrange_basis_poly(data, i)
        result += data.f(x_i) * l * l * (1 - compute_c(data, i) * (x - x_i))

    for i, x_i in enumerate(data.points):
        l = create_lagrange_basis_poly(data, i)
        result += data.f_prime(x_i) * l * l * (x - x_i)

    return simplify(result).collect(x)


if __name__ == "__main__":
    f = Lambda(x, ((x - 1) / 3) ** 8)
    data = HermitData([k / 4 for k in range(5)], f)
    print(create_hermit_poly(data))
