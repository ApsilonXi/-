import math
import matplotlib.pyplot as plt

# Константы
L = 100  # длина стержня
alpha = 1  # коэффициент температуропроводности
n_max = 100
dx = 1  # шаг по пространству
dt = 0.1  # шаг по времени
N = int(L / dx)  # количество пространственных шагов
T_max = 100  # максимальное время
x_values = [i * dx for i in range(N + 1)]

# Задаем начальную функцию
def initial_distribution(x):
    if 50 <= x <= 60:
        return -0.01 * (x - 60)**2 + 1
    elif 60 < x <= 70:
        return (70 - x) / 10
    return 0

# Функция для численного интегрирования методом прямоугольников
def integrate(f, a, b, steps=1000):
    step_size = (b - a) / steps
    total = 0
    for i in range(steps):
        x = a + i * step_size
        total += f(x) * step_size
    return total

# Вычисление коэффициентов Фурье
def compute_fourier_coefficients(n_max):
    coefficients = []
    for m in range(1, n_max + 1):
        # Предварительные вычисления значений углов
        pi = math.pi
        cos_mpi_60 = math.cos(m * pi * 60 / L)
        sin_mpi_60 = math.sin(m * pi * 60 / L)
        cos_mpi_50 = math.cos(m * pi * 50 / L)
        cos_mpi_70 = math.cos(m * pi * 70 / L)
        sin_mpi_70 = math.sin(m * pi * 70 / L)
        
        # Общие члены для сокращения кода
        sin_10m_pi = math.sin(10 * m * pi / L)
        cos_10m_pi = math.cos(10 * m * pi / L)

        # Интеграл I1
        I1_cos_part = -0.01 * cos_mpi_60 * (
            (20 * L**2 * cos_10m_pi) / (m**2 * pi**2) 
            - (2 * L**3 * sin_10m_pi) / (m**3 * pi**3) 
            - (100 * L * sin_10m_pi) / (m * pi)
        )
        I1_sin_part = -0.01 * sin_mpi_60 * (
            (-20 * L**2 * sin_10m_pi) / (m**2 * pi**2) 
            + (2 * L**3 * cos_10m_pi) / (m**3 * pi**3) 
            - (100 * L * cos_10m_pi) / (m * pi)
        )
        I1_cos_50_60 = -(L / (m * pi)) * (cos_mpi_60 - cos_mpi_50)

        # Итоговый I1
        I1 = I1_cos_part + I1_sin_part + I1_cos_50_60

        # Интеграл I2
        I2_cos_part = -(L / (10 * m * pi)) * (70 * cos_mpi_70 - 60 * cos_mpi_60)
        I2_sin_part = (L**2 / (10 * m**2 * pi**2)) * (sin_mpi_70 - sin_mpi_60)

        # Итоговый I2
        I2 = I2_cos_part + I2_sin_part

        # Общий результат
        coefficients.append((2/L) * (I1 + I2))
    return coefficients

b_coefficients = compute_fourier_coefficients(n_max)

# Функция для аналитического решения через ряд Фурье
def u_fourier(x, t, n_max, b_coefficients):
    sum_ = 0
    for n in range(1, n_max + 1):
        bn = b_coefficients[n - 1]
        sum_ += bn * math.exp(-alpha * (n * math.pi / L)**2 * t) * math.sin(n * math.pi * x / L)
    return sum_

# Явное решение задачи
def solve_explicit_method(dx, dt, alpha, N, T_max):
    u = [[0] * (N + 1) for _ in range(int(T_max / dt) + 1)]
    # Начальное условие
    for i in range(N + 1):
        u[0][i] = initial_distribution(i * dx)
    
    # Коэффициент схемы
    r = alpha * dt / dx**2
    
    # Явная схема конечных разностей
    for n in range(0, int(T_max / dt)):
        for i in range(1, N):
            u[n + 1][i] = u[n][i] + r * (u[n][i + 1] - 2 * u[n][i] + u[n][i - 1])
    
    return u

# Решение методом прогонки для неявной схемы
def solve_implicit_method(dx, dt, alpha, N, T_max):
    u = [[0] * (N + 1) for _ in range(int(T_max / dt) + 1)]
    # Начальное условие
    for i in range(N + 1):
        u[0][i] = initial_distribution(i * dx)
    
    r = alpha * dt / dx**2
    
    # Неявная схема через метод прогонки
    for n in range(0, int(T_max / dt)):
        # Правая часть (вектор d)
        d = u[n][1:N]  # Мы не учитываем границы
        
        # Прямой ход прогонки
        P = [0] * (N - 1)
        Q = [0] * (N - 1)
        
        # Параметры для системы уравнений
        A = [-r] * (N - 2)
        B = [(1 + 2 * r)] * (N - 1)
        C = [-r] * (N - 2)

        # Прямой ход
        P[0] = C[0] / B[0]
        Q[0] = d[0] / B[0]
        
        for i in range(1, N - 2):
            denominator = B[i] - A[i - 1] * P[i - 1]
            P[i] = C[i] / denominator
            Q[i] = (d[i] - A[i - 1] * Q[i - 1]) / denominator
        
        # Обратный ход
        u[n + 1][1] = Q[-1]  # Первая неизвестная
        for i in range(N - 3, -1, -1):
            u[n + 1][i + 1] = Q[i] - P[i] * u[n + 1][i + 2]
        
        # Граничные условия
        u[n + 1][0] = u[n + 1][1]  # Левое условие (например, постоянная температура)
        u[n + 1][-1] = u[n + 1][-2]  # Правое условие (например, диффузия)

    return u

# Решаем задачу явным методом
u_explicit = solve_explicit_method(dx, dt, alpha, N, T_max)

# Решаем задачу неявным методом
u_implicit = solve_implicit_method(dx, dt, alpha, N, T_max)

# Моменты времени для вывода
times = [0, 10, 25, 50, 75, 100]
time_steps = [int(t / dt) for t in times]

# Создаем 3 графика на одной фигуре
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 8))

# --- Левый график для решения через ряд Фурье ---
ax1.set_title("Аналитическое решение")
ax1.set_xlabel('x')
ax1.set_ylabel('u(x, t)')
ax1.grid(True)

# График начальной функции
ax1.plot(x_values, [initial_distribution(x) for x in x_values], label='Начальная', linestyle='-', color='black')

# Графики для решения через ряд Фурье с пунктиром
for t in times[1:]:
    u_fourier_values = [u_fourier(x, t, n_max, b_coefficients) for x in x_values]
    ax1.plot(x_values, u_fourier_values, label=f't = {t}', linestyle='-.')

ax1.legend()

# --- Средний график для явного решения ---
ax2.set_title("Явное решение")
ax2.set_xlabel('x')
ax2.set_ylabel('u(x, t)')
ax2.grid(True)

# График начальной функции для явного метода
ax2.plot(x_values, [initial_distribution(x) for x in x_values], label='Начальная', linestyle='-', color='black')

# Графики для явного метода с пунктиром
for idx, t in enumerate(time_steps[1:]):
    ax2.plot(x_values, u_explicit[t], label=f't = {times[idx+1]}', linestyle='-.')

ax2.legend()

# --- Правый график для неявного решения ---
ax3.set_title("Неявное решение")
ax3.set_xlabel('x')
ax3.set_ylabel('u(x, t)')
ax3.grid(True)

# График начальной функции для неявного метода
ax3.plot(x_values, [initial_distribution(x) for x in x_values], label='Начальная', linestyle='-', color='black')

# Графики для неявного метода с пунктиром
for idx, t in enumerate(time_steps[1:]):
    ax3.plot(x_values, u_implicit[t], label=f't = {times[idx+1]}', linestyle='-.')

ax3.legend()

# Показываем все графики
plt.tight_layout()
plt.show()
