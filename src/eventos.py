"""Localización del evento de estado: cruce de rho_pred con theta.

rho_pred = a*U/(c*mu) cruza theta  <=>  U cruza U* = theta*c*mu/a.
Se localiza el instante dentro del paso (no al final) de dos formas:
  1. bisección propia sobre el paso RK4: se busca s en (0, h] con RK4(U, s) = U*;
  2. solve_ivp con `events` (RK45, rtol = atol = 1e-10).
"""
import numpy as np
from scipy.integrate import solve_ivp

from integradores import paso_rk4


def cruce_biseccion(f, U, h, Uumbral, tol=1e-12):
    """Instante s en (0, h] con RK4(U, s) = Uumbral (bisección). Requiere
    U < Uumbral <= RK4(U, h) (cruce hacia arriba dentro del paso)."""
    lo, hi = 0.0, h
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if paso_rk4(f, U, mid) >= Uumbral:
            hi = mid
        else:
            lo = mid
    return hi


def avanzar_con_evento_rk4(f, U, t0, t1, h, Uumbral):
    """Avanza de t0 a t1 con RK4 (pasos <= h) buscando el cruce hacia arriba de
    Uumbral. Devuelve (t_fin, U_fin, cruzo, (ts, Us)) donde (ts, Us) son los nodos."""
    import math
    n = max(1, math.ceil((t1 - t0) / h - 1e-12))
    hh = (t1 - t0) / n
    ts, Us = [t0], [U]
    t = t0
    for i in range(n):
        Un = paso_rk4(f, U, hh)
        if U < Uumbral <= Un:
            s = cruce_biseccion(f, U, hh, Uumbral)
            ts.append(t + s)
            Us.append(Uumbral)
            return t + s, Uumbral, True, (ts, Us)
        t = t0 + (i + 1) * hh
        U = Un
        ts.append(t)
        Us.append(U)
    return t1, U, False, (ts, Us)


def avanzar_con_evento_ivp(f, U, t0, t1, Uumbral, rtol=1e-10, atol=1e-10):
    """Igual que la anterior pero con solve_ivp(events=...), cruce ascendente."""
    def ev(t, y):
        return y[0] - Uumbral
    ev.terminal = True
    ev.direction = 1
    sol = solve_ivp(lambda t, y: [f(y[0])], (t0, t1), [U], events=ev,
                    rtol=rtol, atol=atol, dense_output=False)
    if sol.t_events[0].size:
        te = float(sol.t_events[0][0])
        return te, Uumbral, True, (list(sol.t), list(sol.y[0]))
    return t1, float(sol.y[0, -1]), False, (list(sol.t), list(sol.y[0]))
