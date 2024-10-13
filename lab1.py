def f(t, y):
    """Функция правой части дифференциального уравнения."""
    return -(y**2) / (1 + t**2)

def explicit_method(t0, y0, tf, h):
    """Явный метод Эйлера."""
    n = int((tf - t0) / h)
    t = t0
    y = y0
    results = [(t, y)]
    for _ in range(n):
        y += h * f(t, y)
        t += h
        results.append((t, y))
    return results

def implicit_method(t0, y0, tf, h):
    """Неявный метод Эйлера."""
    n = int((tf - t0) / h)
    t = t0
    y = y0
    results = [(t, y)]
    for _ in range(n):
        y_new = y
        for _ in range(10):  
            y_new = y + h * f(t + h, y_new)
        y = y_new
        t += h
        results.append((t, y))
    return results

def weighted_method(t0, y0, tf, h):
    """Весовой метод."""
    n = int((tf - t0) / h)
    t = t0
    y = y0
    results = [(t, y)]
    for _ in range(n):
        y_explicit = y + h * f(t, y)
        y_weighted = y + 0.5 * h * (f(t, y) + f(t + h, y_explicit))
        y = y_weighted
        t += h
        results.append((t, y))
    return results

def analytical_solution(t, y0):
    """Аналитическое решение для y(t)"""
    # Решение полученное через разделение переменных
    # Уравнение: 1/y = C + arctan(t)
    # Здесь C = 1/y0 - arctan(t0) (если t0 = 0)
    from math import atan
    C = 1 / y0 - atan(0)
    return 1 / (C + atan(t))

# Основной код для тестирования методов
if __name__ == "__main__":
    t0 = 0.0        # Начальное значение t
    y0 = 1.0        # Начальное значение y
    tf = 5.0        # Конечное значение t
    h = 0.1         # Шаг

    # Численные решения
    explicit_results = explicit_method(t0, y0, tf, h)
    implicit_results = implicit_method(t0, y0, tf, h)
    weighted_results = weighted_method(t0, y0, tf, h)

    # Аналитическое решение
    analytical_results = [(t, analytical_solution(t, y0)) for t in [t0 + i * h for i in range(int((tf - t0) / h) + 1)]]

    print("Явный метод Эйлера:")
    for t, y in explicit_results:
        print(f"t = {t:.2f}, y = {y:.5f}")

    print("\nНеявный метод Эйлера:")
    for t, y in implicit_results:
        print(f"t = {t:.2f}, y = {y:.5f}")

    print("\nВесовой метод:")
    for t, y in weighted_results:
        print(f"t = {t:.2f}, y = {y:.5f}")

    print("\nАналитическое решение:")
    for t, y in analytical_results:
        print(f"t = {t:.2f}, y = {y:.5f}")