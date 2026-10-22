"""Capa SD: usuarios activos U(t) (modelo de Bass con abandono).

    dU/dt = (p + q*U/M)(M - U) - delta*U
    delta = delta0 + kW*max(0, W - Wok) + kR*r

Primer término: adopción por innovación (p) y por boca a boca (q*U/M), como el
modelo de Bass de la semana 3. Segundo: abandono, que crece con la latencia W y
con la fracción rechazada r.
"""
from dataclasses import dataclass, replace


@dataclass(frozen=True)
class Params:
    # Capa SD (tiempo en días)
    p: float = 0.002
    q: float = 0.20            # grupo 5: q = 0.20 + 0.05*((5-1) mod 4)
    M: float = 20000.0
    U0: float = 100.0
    delta0: float = 0.003
    kW: float = 1.0            # por día por segundo
    kR: float = 0.5            # por día
    Wok: float = 0.75          # s
    # Capa DES (tiempo en segundos)
    a: float = 0.0005          # solicitudes/s por usuario  (lambda = a*U)
    mu: float = 2.0            # 1/0.5 s
    K: int = 60
    T_des: float = 720.0
    warm: float = 120.0
    # Autoescalado
    theta: float = 0.7
    c0: int = 4
    k_pedido: int = 2          # servidores por pedido
    D: float = 3.0             # días de entrega (grupo 5: 3 + ((5-1) mod 4))
    T: float = 100.0           # horizonte (días)

    def con(self, **kw):
        return replace(self, **kw)


def delta(W, r, P):
    """Tasa de abandono delta = delta0 + kW*max(0, W-Wok) + kR*r."""
    return P.delta0 + P.kW * max(0.0, W - P.Wok) + P.kR * r


def rhs(U, d, P):
    """Lado derecho de la EDO con delta = d constante: dU/dt."""
    return (P.p + P.q * U / P.M) * (P.M - U) - d * U


def umbral_U(c, P):
    """U en el que rho_pred = a*U/(c*mu) = theta:  U* = theta*c*mu/a."""
    return P.theta * c * P.mu / P.a
