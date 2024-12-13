import math
import matplotlib.pyplot as plt
import numpy as np

# Задаем параметры задачи
t = 0.5  # Временной шаг
h = 1    # Пространственный шаг
k = 1    # Коэффициент теплопроводности
x0 = 0   # Левая граница графика
xl = 100 # Правая граница графика
interval = 100 # Количество временных шагов

def initial_distribution(x):
    if 50 <= x <= 60:
        return -0.01 * (x - 60)**2 + 1
    elif 60 < x <= 70:
        return (70 - x) / 10
    return 0

f = initial_distribution      # Начальная функция

# Явный метод для решения задачи теплопроводности
def explicit(f, t, h, k, n, x0, xl):
    x = np.arange(x0, xl + h, h, dtype=float)
    q_prev = np.array([f(_x) for _x in x])
    size = x.size
    for _ in range(1, n):
        q = np.empty(size)
        for i in range(size):
            q_n_i = q_prev[i]
            q_n_ip1 = q_prev[i + 1] if (i + 1 < size) else 0
            q_n_im1 = q_prev[i - 1] if (i - 1 >= 0) else 0
            q[i] = q_n_i + k * t * (q_n_ip1 - 2 * q_n_i + q_n_im1) / (h * h)
        q_prev = q
    return q

# Неявный метод для решения задачи теплопроводности
def implicit(f, t, h, k, n, x0, xl):
    x = np.arange(x0, xl + h, h, dtype=float)
    q_prev = np.array([f(_x) for _x in x])
    size = x.size
    a = b = k / (h * h)
    c = 1 / t + 2 * k / (h * h)
    for _ in range(1, n):
        q = np.empty(size)
        F = [q_i / t for q_i in q_prev]
        alpha = np.zeros(size)
        betta = np.zeros(size + 1)
        z = np.zeros(size)
        # Прямой ход метода прогонки
        for j in range(size - 1):
            z[j] = (c - a * alpha[j])
            alpha[j + 1] = b / (z[j])
            betta[j + 1] = (F[j] + betta[j] * a) / z[j]
        # Обратный ход метода прогонки
        q[size - 1] = betta[size]
        for i in range(size - 1, 0, -1):
            q[i - 1] = betta[i] + alpha[i] * q[i]
        q_prev = q
    return q

'''# Функция для численного интегрирования методом прямоугольников
def integrate(f, a, b, steps=1000):
    step_size = (b - a) / steps
    total = 0
    for i in range(steps):
        x = a + i * step_size
        total += f(x) * step_size
    return total'''

# для аналитического решения задачи
def fourier(f, t, h, k, n, x0, xl):
    x = np.arange(x0, xl + h, h, dtype=float)
    n *= t
    L = xl - x0
    W = math.pi / L
    size = x.size
    C = np.empty(size)

    for m in range(size):
        if m == 0:
            C[m] = 0  
        else:
            integral1 = ((-0.01 ) * (60 - 50) ** 3 + (60 - 50) +1) * (math.sin(W * m * 60) - math.sin(W * m * 50))
            integral2 = (1 / 10) * ((70 - 60) * math.sin(W * m * 70) - (70 - 60) * math.sin(W * m * 60))  

            C[m] = (2 / L) * (integral1 + integral2)

    '''coefficients = []
    for m in range(1, n_max + 1):
        pi = math.pi
        cos_mpi_60 = math.cos(m * pi * 60 / L)
        sin_mpi_60 = math.sin(m * pi * 60 / L)
        cos_mpi_50 = math.cos(m * pi * 50 / L)
        cos_mpi_70 = math.cos(m * pi * 70 / L)
        sin_mpi_70 = math.sin(m * pi * 70 / L)
        
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
        coefficients.append((2/L) * (I1 + I2))'''

    for m in range(1, size):
        integral = 0
        for i in range(len(x) - 1):
            integral += (f(x[i]) * math.sin(W * m * x[i]) + f(x[i + 1]) * math.sin(W * m * x[i + 1])) * (x[i + 1] - x[i]) / 2
        C[m] += (2 / L) * integral

    
    
    q = np.empty(size)
    for i in range(size):
        q[i] = sum(C[m] * math.e ** (-k * W ** 2 * m ** 2 * n) * math.sin(W * m * x[i]) for m in range(1, size))
    return q


x = np.arange(x0, xl + h, h, dtype=float)
initial = np.array([f(_x) for _x in x])

schemas = {
    'Аналитическое': {'func': fourier, 'style': 'dotted', 'color': 'blue'},
    'Явная': {'func': explicit, 'style': 'dashdot', 'color': 'red'},
    'Прогонка': {'func': implicit, 'style': 'dashed', 'color': 'green'}
}

fig, axs = plt.subplots(3, 1, figsize=(7, 14))
for idx, (name, config) in enumerate(schemas.items()):
    func = config['func']   # Получаем функцию
    style = config['style'] # Стиль линии
    color = config['color'] # Цвет линии
    result = func(f, t, h, k, interval, x0, xl)
    axs[idx].plot(x, initial, 'purple', label="Задано")
    axs[idx].plot(x, result, linestyle=style, color=color, label=name)
    axs[idx].set_xlim(0, 100)
    axs[idx].set_ylim(-0.1, 1.1)
    axs[idx].grid(True)
    axs[idx].legend()
    axs[idx].set_title(f"{name}")

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()
