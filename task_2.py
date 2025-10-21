from sympy import *

x = symbols('x')

class NewtonData:
    __slots__ = ('f','n', 'h', 'a')

    def __init__(self, f, n: int, h: float, a: float):
        self.f = f
        self.n = n
        self.h = h
        self.a = a

    def point(self, k: int):
        return self.a + self.h * k

def finite_difference(data: NewtonData, m: int, k: int):
    result = 0

    for j in range(m + 1):
        result += (-1)**(m-j) * binomial(m, j) * data.f(data.point(j))

    return result

def omega_newton(data: NewtonData, k: int):
    result = 1

    for i in range(k):
        result *= (x - data.a)/data.h - i
    
    return result

def create_newton_poly(data: NewtonData):
    result = 0

    for k in range(data.n):
        result += finite_difference(data, k, 0)/factorial(k) * omega_newton(data, k)

    return simplify(result).collect(x)

if __name__ == '__main__':
    f = Lambda(x, ((x-1)/3)**8)
    data = NewtonData(f, 5, 1/4, 0)
    r = create_newton_poly(data)
    print(r)
