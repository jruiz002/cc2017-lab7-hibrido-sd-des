"""Capa DES: cola M/M/c/K propia y fórmulas teóricas de M/M/c.

Fórmulas de clase que se usan:
  * Carga ofrecida          A = lambda/mu          (erlangs)
  * Utilización             rho = lambda/(c*mu)    (estable si rho < 1)
  * Erlang C                C(c,A) = P(esperar)
  * Espera en cola          Wq = C/(c*mu - lambda)
  * Latencia                W = Wq + 1/mu
  * Ley de Little           L = lambda_ef * W,  Lq = lambda_ef * Wq,
                            servidores ocupados = lambda_ef / mu
"""
import heapq
from dataclasses import dataclass
from math import factorial

import numpy as np


# --------------------------------------------------------------------------
# Teoría M/M/c
# --------------------------------------------------------------------------
def erlang_c(c, A):
    """Erlang C: C(c,A) = (A^c/c!)(c/(c-A)) / [sum_{k<c} A^k/k! + (A^c/c!)(c/(c-A))].

    Probabilidad de que una llegada tenga que esperar. Requiere A < c.
    """
    if A >= c:
        return 1.0
    suma = sum(A**k / factorial(k) for k in range(c))
    cola = (A**c / factorial(c)) * (c / (c - A))
    return cola / (suma + cola)


def mmc_teoria(lam, c, mu):
    """Métricas teóricas de M/M/c en estado estable (requiere lam < c*mu)."""
    if lam >= c * mu:
        raise ValueError("M/M/c inestable: lam >= c*mu")
    A = lam / mu
    C = erlang_c(c, A)
    Wq = C / (c * mu - lam)          # espera en cola
    W = Wq + 1.0 / mu                # latencia = Wq + 1/mu
    return {
        "A": A, "rho": lam / (c * mu), "C": C,
        "Wq": Wq, "W": W,
        "Lq": lam * Wq, "L": lam * W,          # Little
        "ocupados": lam / mu,                   # Little aplicada a los servidores
    }


def mmck_bloqueo(lam, c, mu, K):
    """Probabilidad de bloqueo exacta de M/M/c/K (proceso de nacimiento y muerte).

    p_n ∝ A^n/n!            para n <= c
    p_n ∝ A^c/c! * rho^(n-c) para c < n <= K,    P_bloqueo = p_K
    """
    A = lam / mu
    rho = lam / (c * mu)
    p = [1.0]
    for n in range(1, K + 1):
        if n <= c:
            p.append(p[-1] * A / n)
        else:
            p.append(p[-1] * rho)
    p = np.array(p)
    p /= p.sum()
    return p[K]


# --------------------------------------------------------------------------
# Simulación de eventos discretos M/M/c/K
# --------------------------------------------------------------------------
@dataclass
class ResultadoDES:
    """Resumen de una réplica del DES."""
    W: float            # latencia media (aceptadas, llegadas en la parte medida)
    Wq: float           # espera media en cola
    r: float            # fracción rechazada (llegadas medidas)
    util: float         # utilización media por servidor en la parte medida
    L: float            # promedio temporal de clientes en el sistema
    Lq: float           # promedio temporal de clientes en cola
    ocupados: float     # promedio temporal de servidores ocupados
    lam_ef: float       # tasa de aceptación en la parte medida (1/s)
    n_llegadas: int     # llegadas totales de la ventana completa
    n_atendidas: int    # clientes que terminaron servicio dentro de la ventana
    n_rechazadas: int   # llegadas rechazadas (ventana completa)
    n_dentro: int       # clientes que siguen en el sistema al final
    n_medidas: int      # llegadas en la parte medida (t >= calentamiento)


def simular_mmck(lam, c, mu, K, T=720.0, warm=120.0, rng=None, seed=None):
    """Simula M/M/c/K con FIFO durante T segundos y mide después de `warm`.

    Método (eventos discretos con lista de eventos implícita):
      * las llegadas son un proceso de Poisson: tiempos entre llegadas Exp(lam);
      * cada servidor libre se atiende en orden de llegada (FIFO): el inicio de
        servicio es max(llegada, instante en que se libera el primer servidor);
      * si hay K clientes en el sistema al llegar, la solicitud se rechaza.
    Se cuenta cada cliente por su instante de llegada: las métricas de la parte
    medida usan solo las llegadas con t >= warm (calentamiento descartado).
    Los promedios temporales L, Lq y ocupados se integran en [warm, T].
    """
    if rng is None:
        rng = np.random.default_rng(seed)

    servidores = [0.0] * c            # instante en que se libera cada servidor
    heapq.heapify(servidores)
    en_sistema = []                   # heap con instantes de salida de los aceptados

    llegadas = rechazadas = atendidas = 0
    med_llegadas = med_rechazadas = 0
    suma_W = suma_Wq = 0.0
    n_acept_med = 0
    area_L = area_Lq = area_B = 0.0   # integrales sobre [warm, T]

    def solapa(a, b):
        """Longitud de [a,b] ∩ [warm,T]."""
        return max(0.0, min(b, T) - max(a, warm))

    t = 0.0
    while True:
        t += rng.exponential(1.0 / lam)
        if t > T:
            break
        llegadas += 1
        # salidas ocurridas hasta t
        while en_sistema and en_sistema[0] <= t:
            heapq.heappop(en_sistema)
            atendidas += 1
        medida = t >= warm
        if medida:
            med_llegadas += 1
        if len(en_sistema) >= K:                 # sistema lleno: se rechaza
            rechazadas += 1
            if medida:
                med_rechazadas += 1
            continue
        s = rng.exponential(1.0 / mu)
        inicio = max(t, heapq.heappop(servidores))
        salida = inicio + s
        heapq.heappush(servidores, salida)
        heapq.heappush(en_sistema, salida)
        # tiempos del cliente
        area_L += solapa(t, salida)
        area_Lq += solapa(t, inicio)
        area_B += solapa(inicio, salida)
        if medida:
            n_acept_med += 1
            suma_W += salida - t
            suma_Wq += inicio - t

    # Cierre de la ventana: contar las salidas que faltan hasta T
    while en_sistema and en_sistema[0] <= T:
        heapq.heappop(en_sistema)
        atendidas += 1
    dentro = len(en_sistema)

    dur = T - warm
    nan = float("nan")
    return ResultadoDES(
        W=suma_W / n_acept_med if n_acept_med else nan,
        Wq=suma_Wq / n_acept_med if n_acept_med else nan,
        r=med_rechazadas / med_llegadas if med_llegadas else 0.0,
        util=area_B / (c * dur),
        L=area_L / dur, Lq=area_Lq / dur, ocupados=area_B / dur,
        lam_ef=n_acept_med / dur,
        n_llegadas=llegadas, n_atendidas=atendidas,
        n_rechazadas=rechazadas, n_dentro=dentro, n_medidas=med_llegadas,
    )


def conteo_llegadas(lam, T=720.0, warm=120.0, rng=None):
    """Número de llegadas de Poisson(lam) en la parte medida [warm, T]."""
    n_total = 0
    t = 0.0
    n = 0
    while True:
        t += rng.exponential(1.0 / lam)
        if t > T:
            break
        if t >= warm:
            n += 1
    return n
