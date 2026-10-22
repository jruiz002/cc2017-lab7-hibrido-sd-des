"""Integradores propios para dU/dt = f(U) (EDO escalar autónoma).

  * Euler explícito: U_{n+1} = U_n + h f(U_n)                    (orden 1)
  * Runge-Kutta 4:   U_{n+1} = U_n + h/6 (k1 + 2k2 + 2k3 + k4)   (orden 4)
    k1 = f(U), k2 = f(U + h k1/2), k3 = f(U + h k2/2), k4 = f(U + h k3)
"""
import math

import numpy as np


def paso_euler(f, U, h):
    return U + h * f(U)


def paso_rk4(f, U, h):
    k1 = f(U)
    k2 = f(U + 0.5 * h * k1)
    k3 = f(U + 0.5 * h * k2)
    k4 = f(U + h * k3)
    return U + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def integrar(f, U0, t0, t1, h, metodo="rk4"):
    """Integra de t0 a t1 con pasos iguales de tamaño <= h (último paso incluido).

    Devuelve (t, U) en los nodos de la malla. n = ceil((t1-t0)/h) pasos de
    tamaño (t1-t0)/n, así la malla siempre termina exactamente en t1.
    """
    paso = paso_rk4 if metodo == "rk4" else paso_euler
    n = max(1, math.ceil((t1 - t0) / h - 1e-12))
    hh = (t1 - t0) / n
    t = np.linspace(t0, t1, n + 1)
    U = np.empty(n + 1)
    U[0] = U0
    for i in range(n):
        U[i + 1] = paso(f, U[i], hh)
    return t, U
