import numpy as np
import matplotlib.pyplot as plt

# Параметры задачи
dx = 1.0  # Шаг по пространству
dt = 0.1  # Шаг по времени
N = 120  # количество узлов сетки
T = 100  # Время моделирования
u = -0.5  # Скорость конвекции
mu = 1.0  # Коэффициент диффузии mu
sigma = 0.5  # Параметр схемы, отвечающий за стабильность разностной аппроксимации (обеспечивает устойчивость и точность)
omega = np.pi / 20  # Частота для синусоидальной функции на левой границе
C0 = -2  # Концентрация на правой границе
alpha = 0.1  # Коэффициент аппроксимации правой границы

# Инициализация массивов
q = np.zeros(N + 1)
A = np.zeros(N + 1)  # Будут использоваться для аппроксимации уравнения
B = np.zeros(N + 1)
C = np.zeros(N + 1)
F = np.zeros(N + 1)

for i in range(60, 80):
    q[i] = -(1 / 100) * (i - 70) / 2 + 1

# Метод прогонки
def run_through_method(A, B, C, F):
    zn = np.zeros(len(F))
    alf = np.zeros(len(F))
    bet = np.zeros(len(F))
    q = np.zeros(len(F))

    # Прямой ход (вычисление вспомогательных массивов alf и bet для оптимального решения)
    zn[0] = C[0] if C[0] != 0 else 1e-10
    alf[1] = B[0] / zn[0]
    bet[1] = F[0] / zn[0]

    for i in range(1, N):
        zn[i] = C[i] - alf[i] * A[i]
        zn[i] = zn[i] if zn[i] != 0 else 1e-10
        alf[i + 1] = B[i] / zn[i]
        bet[i + 1] = (F[i] + A[i] * bet[i]) / zn[i]

    q[N] = (F[N] + A[N] * bet[N - 1]) / (C[N] - A[N] * alf[N - 1])

    # Обратный ход (построение решения q на основе вспомогательных массивов)
    for i in range(N - 1, -1, -1):
        q[i] = alf[i + 1] * q[i + 1] + bet[i + 1]

    return q

# Временной цикл
time_cycle = np.arange(0, T + dt, dt)
q_results = []
for t in time_cycle:
    # левое: синусоидальная зависимость от времени
    q[0] = np.sin(omega * t)

    # для внутренних точек
    for i in range(1, N):
        A[i] = (u / (2 * dx) + mu / (dx * dx)) * sigma
        B[i] = (-u / (2 * dx) + mu / (dx * dx)) * sigma
        C[i] = 1.0 / dt + A[i] + B[i] #д

    A[N] = (u / (2 * dx) + mu / (dx * dx)) * sigma
    B[N] = 0
    C[N] = 1.0 / (2 * dt) + B[N]

    # для внутренних узлов
    for i in range(1, N):
        F[i] = q[i] * (1 / dt - 2 * (1 - sigma) * mu / (dx * dx)) + q[i + 1] * (-u / (2 * dx) + mu / (dx * dx)) * (1 - sigma) + q[i - 1] * (u / (2 * dx) + mu / (dx * dx)) * (1 - sigma)
    # для последнего узла
    F[N] = q[N] * (1 / (2 * dt) - (-u / (2 * dx) + mu / (dx * dx)) * (1 - sigma)) + q[N - 1] * (-u / (2 * dx) + mu / (dx * dx)) * (1 - sigma)

    q = run_through_method(A, B, C, F)

    q_results.append(q.copy())

q_results = {int(t): q for t, q in zip(time_cycle, q_results)}

# Визуализация результатов
plt.figure(figsize=(10, 6))
x_coords = np.arange(0, N + 1)
initial_function = np.maximum(-((x_coords - 70)  2) / 100 + 1, 0)

# График исходной функции
plt.plot(x_coords, initial_function, 'k--', label="Парабола")

# Добавление графиков для каждого момента времени
colors = ["#FF6F61", "#6B5B95", "#88B04B", "#F7CAC9", "#92A8D1", "#955251"]
time_points = [0, 10, 25, 50, 75, 100]

for countT, t in enumerate(time_points):
    plt.plot(x_coords, q_results[t], label=f"t = {t} c", color=colors[countT % len(colors)])

plt.xlabel("x")
plt.ylabel("y")
plt.xlim(0, 100)
plt.ylim(-0.5, 1.5)
plt.title("Решение задачи диффузии-конвекции с граничным условием")
plt.legend()
plt.grid()

# Сохранение графика
plt.savefig("diffusion_conveection_solution.png", dpi=300, bbox_inches='tight') # e
plt.show()