from sympy import *
from sympy.abc import x


class SplineData:
    slots = ("f", "xs")

    def __init__(self, f, xs):
        self.f = f
        self.xs = xs
        # self.n = len(xs)


def create_first_order_base_spline(data: SplineData, k: int):
    xs = data.xs
    # cond3 = (0, x < xs[0])
    # cond4 = (0, x > xs[len(xs)])

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
    cond3 = ((xs[k + 1] - x) / (xs[k + 1] - xs[k]), x <= xs[k + 1])
    cond4 = (0, True)
    return Piecewise(cond1, cond2, cond3, cond4)


def create_first_order_spline(data: SplineData):
    result = 0

    for k, x_k in enumerate(data.xs):
        b = create_first_order_base_spline(data, k)
        print(b)
        result += data.f(x_k) * b

    # return piecewise_exclusive(piecewise_fold(result))
    return result


def equispaced_nodes(a, b, n):
    return [a + i * (b - a) / (n - 1) for i in range(n)]


def chebyshev_nodes(a, b, n):
    return [
        (a + b) / 2 + (b - a) * cos((2 * k - 1) * pi / (2 * n)) / 2
        for k in range(1, n + 1)
    ]


def create_third_order_spline(data: SplineData):
    xs = data.xs
    syms = [list(symbols(f"a{k}(0:4)")) for k in range(1, len(xs))]
    # print(syms)
    conds = [(0, x < xs[0])]
    pieces = []

    for i in range(len(xs) - 1):
        p = syms[i][0] + syms[i][1] * x + syms[i][2] * x**2 + syms[i][3] * x**3
        pieces.append(p)

        conds.append((p, x <= xs[i + 1]))

    conds.append((0, True))
    # print(conds)
    # print(pieces)
    # print(pieces[0].diff(x))

    p = Piecewise(*conds)
    # print(p)

    equations = []

    for i in range(len(xs)):
        equations.append(Eq(p.subs(x, xs[i]), data.f(xs[i])))

    for i in range(len(pieces) - 1):
        equations.append(
            Eq(pieces[i].subs(x, xs[i + 1]), pieces[i + 1].subs(x, xs[i + 1]))
        )
        equations.append(
            Eq(
                pieces[i].diff(x).subs(x, xs[i + 1]),
                pieces[i + 1].diff(x).subs(x, xs[i + 1]),
            )
        )
        equations.append(
            Eq(
                pieces[i].diff(x, 2).subs(x, xs[i + 1]),
                pieces[i + 1].diff(x, 2).subs(x, xs[i + 1]),
            )
        )

    equations.append(Eq(pieces[0].diff(x, 2).subs(x, xs[0]), 0))
    equations.append(Eq(pieces[-1].diff(x, 2).subs(x, xs[-1]), 0))

    for k, v in solve(equations).items():
        p = p.subs(k, v)

    return p


if __name__ == "__main__":
    a, b, n = -1, 1, 9
    data = SplineData(Lambda(x, sin(x)), equispaced_nodes(a, b, n))
    # f = create_first_order_spline(data)
    # plot(f, (x, a - 1, b + 1), show=True)
    f = create_third_order_spline(data)
    plot(f, (x, a - 0.5, b + 0.5), show=True)
    pass
