import math
import matplotlib.pyplot as plt

def f(t, y):
    return -(y**2)/math.sqrt(4+9*t**2)

# Аналитическое решение
def analytical_solution(t):
    return 1/(math.log1p(math.sqrt(1+t**2)+t)+1)

# Явная схема
def explicit_scheme(t0, y0, h, N):
    t = [t0 + i*h for i in range(N + 1)]
    y = [y0]
    for i in range(1, N + 1):
        y.append(y[i - 1] + h * f(t[i - 1], y[i - 1]))
    return t, y

# Неявная схема
def implicit_scheme(t0, y0, h, N):
    t = [t0 + i*h for i in range(N + 1)]
    y = [y0]
    for i in range(1, N + 1):
        y_next = y[i - 1]  
        for _ in range(10):  
            y_next = y[i - 1] + (h / 2) * (f(t[i], y_next) + f(t[i - 1], y[i - 1]))
        y.append(y_next)
    return t, y

# Весовая схема
def weighted_scheme(t0, y0, h, N):
    t = [t0 + i*h for i in range(N + 1)]
    y = [y0]
    for i in range(1, N + 1):
        y.append(y[i - 1] + (h / 2) * (f(t[i - 1], y[i - 1]) + f(t[i], y[i - 1])))
    return t, y

if __name__ == "__main__":
    t0 = 0        
    y0 = 1        
    h = 0.1       # шаг
    N = 100       # количество шагов

    t_analytical = [t0 + i*h for i in range(N + 1)]
    y_analytical = [analytical_solution(t) for t in t_analytical]

    t_explicit, y_explicit = explicit_scheme(t0, y0, h, N)
    t_implicit, y_implicit = implicit_scheme(t0, y0, h, N)
    t_weighted, y_weighted = weighted_scheme(t0, y0, h, N)

    print("t\tAnalytical\tExplicit\tImplicit\tWeighted")
    for i in range(len(t_explicit)):
        print(f"{y_analytical[i]:.4f}\t{t_explicit[i]:.2f}\t\t{y_explicit[i]:.4f}\t\t{y_implicit[i]:.4f}\t\t{y_weighted[i]:.4f}")

    plt.figure(figsize=(10, 6))
    plt.plot(t_explicit, y_explicit, label='Explicit Scheme', marker='o')
    plt.plot(t_implicit, y_implicit, label='Implicit Scheme', marker='x')
    plt.plot(t_weighted, y_weighted, label='Weighted Scheme', marker='s')
    plt.plot(t_analytical, y_analytical, label='Analytical Solution', linestyle='dashed')

    plt.xlabel('t')
    plt.ylabel('y(t)')
    plt.legend()
    plt.grid()
    plt.show()