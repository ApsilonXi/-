import numpy as np
import matplotlib.pyplot as plt

L, T, Nx, Nt, v = 150, 50, 90, 200, -0.5
dx, dt = L / Nx, T / Nt
x, t = np.linspace(0, L, Nx), np.linspace(0, T, Nt)

def wave(x):
    if 70 <= x <= 90:
        return np.sqrt(1 - ((x - 80) / 10) ** 2)
    else:
        return 0

q_initial = np.array([wave(xi) for xi in x])
q_old = q_initial.copy() 
def center():
    q = q_initial.copy()
    for n in range(Nt - 1):
        q_new = q.copy()
        for i in range(1, Nx - 1):
            q_new[i] = q[i] - v * dt / (2 * dx) * (q[i + 1] - q[i - 1]) + 0.2 * (q[i + 1] - 2 * q[i] + q[i - 1])
        q = q_new
    return q
def kabare():
    q = q_initial.copy()
    global q_old
    for n in range(Nt - 1):
        q_new = q.copy()
        for i in range(1, Nx - 1):
            q_new[i] =  q_new[i] = q[i] - q[i + 1] + q_old[i + 1] - v * (2 * dt) / dx * (q[i + 1] - q[i])
        q_old = q.copy()
        q = q_new
    return q

def left():
    q = q_initial.copy()
    for n in range(Nt - 1):
        q_new = q.copy()
        for i in range(Nx - 1):
            q_new[i] = q[i] - v * dt / dx * (q[i + 1] - q[i])
        q = q_new
    return q

def mixed():
    q = q_initial.copy()
    global q_old
    for n in range(Nt - 1):
        q_new = q.copy()
        for i in range(1, Nx - 1):
            if n % 2 == 0:  # центральная разностная схема на четных шагах
                q_new[i] = q[i] - v * dt / (2 * dx) * (q[i + 1] - q[i - 1])
            else:  # схема Кабаре на нечетных шагах
                q_new[i] = q[i] - q[i + 1] + q_old[i + 1] - v * (2 * dt) / dx * (q[i + 1] - q[i])
        q_old = q.copy()
        q = q_new
    return q


def cross_scheme():
    q_old = q_initial.copy()  
    q = q_initial.copy()    
    q_new = np.zeros_like(q)  

    for n in range(1, Nt):    
        for i in range(1, Nx - 1):
            q_new[i] = q_old[i] - v * dt / dx * (q[i + 1] - q[i - 1])
        q_old = q.copy()
        q = q_new.copy()

    return q

q_center = center()
q_kabare = kabare()
q_left = left()
q_smes = mixed()
q_cross = cross_scheme()
fig, axs = plt.subplots(2, 2, figsize=(11, 6))

# центральная схема
axs[0, 0].plot(x, q_center, 'b-', label='Центральная схема', linewidth=2, color='black')
axs[0, 0].set_xlabel('x')
axs[0, 0].set_ylabel('q')
axs[0, 0].grid(True)
axs[0, 0].legend()

# схема Кабаре
axs[0, 1].plot(x, q_kabare, 'g-', label='Схема Кабаре', linewidth=2, color='green')
axs[0, 1].set_xlabel('x')
axs[0, 1].set_ylabel('q')
axs[0, 1].grid(True)
axs[0, 1].legend()

# левый уголок
axs[1, 0].plot(x, q_left, 'm-', label='Левый уголок', linewidth=2, color='grey')
axs[1, 0].set_xlabel('x')
axs[1, 0].set_ylabel('q')
axs[1, 0].grid(True)
axs[1, 0].legend()

# смешанная схема
axs[1, 1].plot(x, q_smes, 'c-', label='Смешанная схема', linewidth=2, color='red')
axs[1, 1].set_xlabel('x')
axs[1, 1].set_ylabel('q')
axs[1, 1].grid(True)
axs[1, 1].legend()

plt.tight_layout()
plt.show()