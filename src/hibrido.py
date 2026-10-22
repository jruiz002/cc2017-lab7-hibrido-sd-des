"""Híbrido SD-DES con autoescalado.

Acoplamiento (cada H días):  SD --(lambda = a*U, c)--> DES --(W, r)--> SD.
W y r se mantienen constantes hasta el siguiente intercambio, de modo que dentro
de cada intervalo el SD es una EDO autónoma con delta constante.

Autoescalado (evento de estado + evento de tiempo):
  * estado: rho_pred = a*U/(c*mu) cruza theta hacia arriba sin pedido pendiente
    -> se piden k_pedido servidores;
  * tiempo: los servidores llegan D días después del pedido; si en ese instante
    rho_pred sigue > theta se ordena otro pedido de inmediato.

Dos formas de detectar el cruce del umbral:
  * modo="evento": instante localizado dentro del paso (bisección con RK4 o
    solve_ivp events); el pedido y la llegada ocurren en sus instantes exactos;
  * modo="ingenuo": se revisa rho_pred solo al final de cada paso dt.
"""
import numpy as np

from des import simular_mmck
from eventos import avanzar_con_evento_ivp, avanzar_con_evento_rk4
from integradores import integrar
from sd import delta, rhs, umbral_U

EPS = 1e-9


def llamar_des(U, c, P, seed, k):
    """DES de la hora pico: lambda = a*U. Semilla propia por intercambio (k)."""
    rng = np.random.default_rng([seed, k])
    return simular_mmck(P.a * U, c, P.mu, P.K, P.T_des, P.warm, rng=rng)


def simular_hibrido(P, H=1.0, modo="evento", dt=0.25, lazo=True, seed=0,
                    h_rk4=0.25, localizador="biseccion"):
    """Simula P.T días del híbrido y devuelve un diccionario de resultados.

    lazo=False reproduce la cascada (kW = kR = 0: el DES no afecta a los usuarios).
    modo="evento": h_rk4 es el paso de RK4 (y se localiza el cruce);
    modo="ingenuo": el paso es dt y el umbral se revisa solo al final del paso.
    """
    paso = h_rk4 if modo == "evento" else dt
    t, U, c = 0.0, P.U0, P.c0
    pend = None                       # instante de llegada del pedido pendiente
    pedidos, llegadas = [], []        # (t_pedido, c_antes) / (t_llegada, c_despues)
    traj_t, traj_U = [0.0], [U]
    c_t, c_v = [0.0], [c]
    ex_t, ex_W, ex_r, ex_U, ex_c = [], [], [], [], []
    rho_max = 0.0

    def pedir(tp):
        nonlocal pend
        pend = tp + P.D
        pedidos.append((tp, c))

    k = 0
    while t < P.T - EPS:
        # --- intercambio SD -> DES -> SD en t = k*H ---
        res = llamar_des(U, c, P, seed, k)
        W, r = res.W, res.r
        ex_t.append(t); ex_W.append(W); ex_r.append(r); ex_U.append(U); ex_c.append(c)
        d = delta(W, r, P) if lazo else P.delta0
        f = lambda u, d=d: rhs(u, d, P)
        t_ex = min((k + 1) * H, P.T)
        k += 1

        if modo == "evento":
            while t < t_ex - EPS:
                if pend is None and U >= umbral_U(c, P):   # reorden inmediato
                    pedir(t)
                t_stop = min(t_ex, pend) if pend is not None else t_ex
                if pend is None:
                    if localizador == "biseccion":
                        tf, Uf, cruzo, (ts, Us) = avanzar_con_evento_rk4(
                            f, U, t, t_stop, paso, umbral_U(c, P))
                    else:
                        tf, Uf, cruzo, (ts, Us) = avanzar_con_evento_ivp(
                            f, U, t, t_stop, umbral_U(c, P))
                    traj_t += list(ts[1:]); traj_U += list(Us[1:])
                    t, U = tf, Uf
                    if cruzo:
                        pedir(t)
                        continue
                else:
                    ts, Us = integrar(f, U, t, t_stop, paso)
                    traj_t += list(ts[1:]); traj_U += list(Us[1:])
                    t, U = t_stop, Us[-1]
                if pend is not None and t >= pend - EPS:     # llegada de servidores
                    c += P.k_pedido
                    llegadas.append((t, c)); c_t.append(t); c_v.append(c)
                    pend = None
        else:                                                # versión ingenua
            while t < t_ex - EPS:
                t_next = min(t + dt, t_ex)
                ts, Us = integrar(f, U, t, t_next, paso)
                traj_t += list(ts[1:]); traj_U += list(Us[1:])
                t, U = t_next, Us[-1]
                if pend is not None and t >= pend - EPS:
                    c += P.k_pedido
                    llegadas.append((t, c)); c_t.append(t); c_v.append(c)
                    pend = None
                # revisión de rho_pred solo al final del paso
                if pend is None and U >= umbral_U(c, P):
                    pedir(t)
        rho_max = max(rho_max, P.a * U / (c * P.mu))

    return {
        "t": np.array(traj_t), "U": np.array(traj_U),
        "c_t": np.array(c_t), "c_v": np.array(c_v),
        "ex_t": np.array(ex_t), "ex_W": np.array(ex_W), "ex_r": np.array(ex_r),
        "ex_U": np.array(ex_U), "ex_c": np.array(ex_c),
        "pedidos": pedidos, "llegadas": llegadas,
        "U_final": U, "c_final": c,
        "W_max": float(np.nanmax(ex_W)), "r_max": float(np.max(ex_r)),
        "t_pedido1": pedidos[0][0] if pedidos else float("nan"),
        "t_llegada1": llegadas[0][0] if llegadas else float("nan"),
    }


def rho_fina(R, P, dt=0.01):
    """rho_pred(t) = a*U(t)/(c(t)*mu) remuestreada en una malla fina (dt = 0.01 días).

    Medida continua de la congestión del SD que no depende de que el DES se consulte
    solo en los intercambios ni del paso del integrador: U se interpola linealmente
    entre los nodos de la trayectoria y c(t) es constante a trozos (el servidor nuevo
    cuenta desde su instante de llegada)."""
    t = np.append(np.arange(0.0, P.T, dt), P.T)
    U = np.interp(t, R["t"], R["U"])
    idx = np.searchsorted(R["c_t"], t, side="right") - 1
    return t, P.a * U / (R["c_v"][idx] * P.mu)


def exposicion(R, P, dt=0.01):
    """Resumen de congestión: (rho_pred máx, días con rho_pred > theta, días con rho_pred >= 1)."""
    t, rho = rho_fina(R, P, dt)
    return rho.max(), dt * np.sum(rho[:-1] > P.theta), dt * np.sum(rho[:-1] >= 1.0)
