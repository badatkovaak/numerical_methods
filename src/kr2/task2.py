from sympy import *
from sympy.abc import x
from src.kr2.common import SplineData


def third_order_spline_by_definition(data: SplineData):
    xs = data.xs
    syms = [list(symbols(f"a{k}(0:4)")) for k in range(1, len(xs))]
    conds = [(0, x < xs[0])]
    pieces = []

    for i in range(len(xs) - 1):
        p = syms[i][0] + syms[i][1] * x + syms[i][2] * x**2 + syms[i][3] * x**3
        pieces.append(p)

        conds.append((p, x <= xs[i + 1]))

    conds.append((0, True))

    p = Piecewise(*conds)

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


def third_order_spline_by_moments(data: SplineData):
    xs = data.xs
    syms = [list(symbols(f"a{k}(0:2)")) for k in range(1, len(xs))]
    moments = list(symbols(f"M:{len(xs)}"))
    conds = [(0, x < xs[0])]
    pieces = []

    for i in range(len(xs) - 1):
        piece = (
            moments[i] * (x - xs[i]) ** 2 / 2
            + ((moments[i] - moments[i + 1]) / (xs[i + 1] - xs[i]))
            * ((x - xs[i]) ** 3)
            / 6
            + syms[i][0] * (x - xs[i])
            + syms[i][1]
        )
        pieces.append(piece)
        conds.append((piece, x <= xs[i + 1]))

    p = Piecewise(*conds)

    print(pieces)
    print(p)
    return p
