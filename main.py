import matplotlib.pyplot as plt
import numpy as np
from sympy import *


def draw_function(f, a, b, epsilon):
    xs = np.linspace(a, b, num=100)
    ys = [calculate_integral(f, a, x, epsilon, False) for x in xs]
    plt.plot(xs, ys)
    plt.show()


def central_rectangle_method(f, n, delta, b):
    h = (b - delta) / n
    sum = 0
    k = 0
    while k < n:
        i = k + 0.5
        x_k = delta + (i * h)
        sum += f.subs(x, x_k)
        k += 1
    return h * sum


def calculate_integral(f, a, b, epsilon, verbose = True):
    f1 = Integral(f, (x, a, d))
    f2 = Integral(f, (x, d, b))

    if verbose:
        print("I1 =", f1)
        print("I2 =", f2)
        print()

    g1 = 1 / sqrt(1 + x ** 2)
    g2 = 1 / sqrt(x)

    g1_max = maximum(g1, x, Interval(a, a + 0.0001))
    if g1_max != g1.subs(x, a):
        g1_max = g1.subs(x, d)

    I1 = g1_max * integrate(g2, (x, a, d))

    epsilon_1 = epsilon / 1.3
    epsilon_2 = epsilon - epsilon_1
    if verbose:
        print("epsilon_1 =", epsilon_1)
        print("epsilon_2 =", epsilon_2)
        print("I1 <=", f"max({g1})", "*", Integral(g2, (x, 0, d)), "=", I1)
        print()

    delta = solve(I1 - epsilon_1, d)[0].evalf()
    if verbose:
        print("delta =", delta)
        print()

    n = 2
    I_n = 0
    I_2n = central_rectangle_method(f, n, delta, b)

    while (abs(I_n - I_2n) / 3) > epsilon_2:
        I_n = I_2n
        n *= 2
        I_2n = central_rectangle_method(f, n, delta, b)

    res = I_2n + epsilon_1
    if verbose:
        print(f"n = {n}: |{I_n} - {I_2n}|/ (2^2 - 1) = {abs(I_n - I_2n) / 3} <= {epsilon_2}")
        print("Result = ", I_2n, "+", epsilon_1, "=", res)
    return res


if __name__ == '__main__':
    init_printing(use_unicode=False, wrap_line=False)
    x = Symbol('x')
    d = Symbol('d')
    a = 0
    b = 5
    epsilon = 0.1
    f = 1/sqrt(x*(1+x**2))

    calculate_integral(f, a, b, epsilon)

    # draw_function(f, a, b, epsilon)





