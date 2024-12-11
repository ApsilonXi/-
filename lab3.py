import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Константы
L = 100  # длина стержня
alpha = 1  # коэффициент температуропроводности
n_max = 100
dx = 1  # шаг по пространству
dt = 0.1  # шаг по времени
N = int(L / dx)  # количество пространственных шагов
T_max = 100  # максимальное время
x_values = np.linspace(0, L, N + 1)

def initial_distribution(x):
    if 0 <= x < 50:
        return 0
    elif 50 <= x < 60:
        return ((x - 50)**2) / 100
    elif 60 <= x < 70:
        return 1 - (x - 60) / 10
    else:
        return 0

def compute_fourier_coefficients(n_max):
    coefficients = []
    for n in range(1, n_max + 1):
        integrand = lambda x: initial_distribution(x) * np.sin(n * np.pi * x / L)
        integral, _ = quad(integrand, 0, L)
        b_n = (2 / L) * integral
        coefficients.append(b_n)
    return np.array(coefficients)


b_coefficients = compute_fourier_coefficients(n_max)


def u_fourier(x, t, n_max, b_coefficients): #расчёт
    sum_ = 0
    for n in range(1, n_max + 1):
        bn = b_coefficients[n - 1]
        sum_ += bn * np.exp(-alpha * (n * np.pi / L)**2 * t) * np.sin(n * np.pi * x / L)
    return sum_

# Явное решение задачи
def solve_explicit_method(dx, dt, alpha, N, T_max):
    u = np.zeros((int(T_max / dt) + 1, N + 1))
    # Начальное условие
    for i in range(N + 1):
        u[0, i] = initial_distribution(i * dx)
    
    # Коэффициент схемы
    r = alpha * dt / dx**2
    
    # Явная схема конечных разностей
    for n in range(0, int(T_max / dt)):
        for i in range(1, N):
            u[n + 1, i] = u[n, i] + r * (u[n, i + 1] - 2 * u[n, i] + u[n, i - 1])
    
    return u

# Решаем задачу явным методом
u_explicit = solve_explicit_method(dx, dt, alpha, N, T_max)

# Решение методом прогонки для неявной схемы
def solve_implicit_method(dx, dt, alpha, N, T_max):
    u = np.zeros((int(T_max / dt) + 1, N + 1))
    # Начальное условие
    for i in range(N + 1):
        u[0, i] = initial_distribution(i * dx)
    
    r = alpha * dt / dx**2
    
    # Неявная схема через метод прогонки
    for n in range(0, int(T_max / dt)):
        # Правая часть (вектор d)
        d = u[n, 1:N]  # Мы не учитываем границы
        
        # Прямой ход прогонки
        P = np.zeros(N - 1)
        Q = np.zeros(N - 1)
        
        # Параметры для системы уравнений
        A = -r * np.ones(N - 2)
        B = (1 + 2 * r) * np.ones(N - 1)
        C = -r * np.ones(N - 2)

        # Прямой ход
        P[0] = C[0] / B[0]
        Q[0] = d[0] / B[0]
        
        for i in range(1, N - 2):
            denominator = B[i] - A[i - 1] * P[i - 1]
            P[i] = C[i] / denominator
            Q[i] = (d[i] - A[i - 1] * Q[i - 1]) / denominator
        
        # Обратный ход
        u[n + 1, 1] = Q[-1]  # Первая неизвестная
        for i in range(N - 3, -1, -1):
            u[n + 1, i + 1] = Q[i] - P[i] * u[n + 1, i + 2]
        
        # Граничные условия
        u[n + 1, 0] = u[n + 1, 1]  # Левое условие (например, постоянная температура)
        u[n + 1, -1] = u[n + 1, -2]  # Правое условие (например, диффузия)

    return u

# Решаем задачу неявным методом
u_implicit = solve_implicit_method(dx, dt, alpha, N, T_max)

# Моменты времени для вывода
times = [0, 10, 25, 50, 75, 100]
time_steps = [int(t / dt) for t in times]

# Функция для вычисления RMSE
def compute_rmse(true_values, approx_values):
    return np.sqrt(np.mean((true_values - approx_values) ** 2))

# Графики разностей
fig_diff, axs = plt.subplots(len(time_steps) - 1, 1, figsize=(10, 5 * (len(time_steps) - 1)))
for idx, t in enumerate(time_steps[1:]):
    u_fourier_values = [u_fourier(x, times[idx + 1], n_max, b_coefficients) for x in x_values]
    difference_explicit_fourier = u_explicit[t, :] - u_fourier_values
    difference_implicit_fourier = u_implicit[t, :] - u_fourier_values
    
    # Вычисление RMSE
    rmse_explicit = compute_rmse(u_fourier_values, u_explicit[t, :])
    rmse_implicit = compute_rmse(u_fourier_values, u_implicit[t, :])
    
    # Построение графиков разностей
    axs[idx].plot(x_values, difference_explicit_fourier, label='Явное - Аналитическое', linestyle='--')
    axs[idx].plot(x_values, difference_implicit_fourier, label='Неявное - Аналитическое', linestyle='-.')
    axs[idx].set_title(f'Разность температур для t = {times[idx + 1]}')
    axs[idx].set_xlabel('x')
    axs[idx].set_ylabel('Температура')
    axs[idx].legend()
    axs[idx].grid(True)
    axs[idx].axhline(0, color='black', linewidth=0.5, linestyle='--')  # Линия y=0
    
    # Вывод RMSE
    print(f'RMSE для t = {times[idx + 1]}: Явное = {rmse_explicit:.4f}, Неявное = {rmse_implicit:.4f}')

# Показываем графики разностей
plt.tight_layout()
plt.show()
