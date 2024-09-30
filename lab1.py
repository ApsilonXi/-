import numpy as np
import matplotlib.pyplot as plt

# Параметры
h = 0.1          # Шаг интегрирования
T = 10           # Конечное время
N = int(T / h)   # Количество шагов
t = np.linspace(0, T, N + 1)  # Время

# Изначальные условия
y0 = 1

# Явная схема (Метод Эйлера)
def explicit_euler(h, y0, t):
    y = np.zeros(N + 1)
    y[0] = y0
    for n in range(N):
        y[n + 1] = y[n] - h * (y[n]**2) / (1 + t[n]**2)
    return y

# Неявная схема (Метод Эйлера)
def implicit_euler(h, y0, t):
    y = np.zeros(N + 1)
    y[0] = y0
    for n in range(N):
        # Нелинейное уравнение, решение методом простых итераций
        y_prev = y[n]  # Начальное приближение
        for _ in range(10): # Итерации для нахождения y[n+1]
            y_next = y_prev + h * (-(y_prev**2) / (1 + (t[n] + h)**2))
            y_prev = y_next
        y[n + 1] = y_next
    return y

# Весовая схема (Метод Рунге-Кутты 2-го порядка)
def runge_kutta(h, y0, t):
    y = np.zeros(N + 1)
    y[0] = y0
    for n in range(N):
        k1 = - (y[n]**2) / (1 + t[n]**2)
        k2 = - ((y[n] + h * k1)**2) / (1 + (t[n] + h)**2)
        y[n + 1] = y[n] + (h / 2) * (k1 + k2)
    return y

# Аналитическое решение
def analytical_solution(t, y0):
    return 1 / (np.arctan(t) + 1 / y0)

# Решение
y_explicit = explicit_euler(h, y0, t)
y_implicit = implicit_euler(h, y0, t)
y_rk = runge_kutta(h, y0, t)
y_analytical = analytical_solution(t, y0)

# Вывод результата в консоль
print("Явная схема (Метод Эйлера):", y_explicit)
print("Неявная схема (Метод Эйлера):", y_implicit)
print("Весовая схема (Метод Рунге-Кутты):", y_rk)
print("Аналитическое решение:", y_analytical)

# Построение графика
plt.figure(figsize=(10, 6))
plt.plot(t, y_explicit, label='Явная схема', linestyle='-', marker='o')
plt.plot(t, y_implicit, label='Неявная схема', linestyle='-', marker='x')
plt.plot(t, y_rk, label='Метод Рунге-Кутты', linestyle='-', marker='d')
plt.plot(t, y_analytical, label='Аналитическое решение', linestyle='--')
plt.title('Сравнение численных и аналитического решений')
plt.xlabel('Время (t)')
plt.ylabel('Решение (y)')
plt.legend()
plt.grid()
plt.show()