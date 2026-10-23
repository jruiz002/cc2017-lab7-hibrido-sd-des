"""
Modelo híbrido System Dynamics (SD) + Simulación de Eventos Discretos (DES)
Curso: Modelación y Simulación

- Capa SD (tiempo en días): dinámica de usuarios activos U(t) (modelo tipo Bass con abandono).
- Capa DES (segundos, hora pico): cola M/M/c/K que devuelve latencia W y fracción rechazada r.
- Autoescalado: pedidos de 2 servidores con retraso D = 3 días.
- Acoplamiento: cada H = 1 día el SD le pasa (lambda, c) al DES y el DES devuelve (W, r).
"""

import heapq
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# 1. PARÁMETROS DEL MODELO
# =============================================================================
# --- Capa SD ---
p = 0.002        # coeficiente de innovación (adopción "externa")
q = 0.20         # coeficiente de imitación (boca a boca)
M = 20000        # mercado potencial
U0 = 100         # usuarios activos iniciales
delta0 = 0.003   # tasa base de abandono (1/día)
kW = 1.0         # sensibilidad del abandono a la latencia excesiva
kR = 0.5         # sensibilidad del abandono a los rechazos
Wok = 0.75       # latencia aceptable (s)

# --- Capa DES ---
a = 0.0005       # solicitudes por segundo por usuario
mu = 2.0         # tasa de servicio por servidor (media 0.5 s)
K = 60           # capacidad máxima del sistema (en servicio + en cola)
T_DES = 720.0    # duración de cada corrida DES (s)
WARMUP = 120.0   # calentamiento descartado (s)

# --- Autoescalado ---
c0 = 4           # servidores iniciales
RHO_UMBRAL = 0.7 # umbral de utilización predicha
DELTA_C = 2      # servidores por pedido
D = 3.0          # retraso de entrega (días)

# --- Simulación / acoplamiento ---
T_FINAL = 100.0  # horizonte (días)
H = 1.0          # intervalo de intercambio SD <-> DES (días)
DT = 0.01        # paso de integración del SD (días); H debe ser múltiplo de DT


# =============================================================================
# 2. CAPA DES: COLA M/M/c/K
# =============================================================================
def simular_des(lam, c, rng):
    """
    Simula una cola M/M/c/K FCFS durante T_DES segundos.

    Como la disciplina es FCFS, cada solicitud aceptada empieza a atenderse en
    cuanto se libera el primer servidor, así que podemos calcular su salida de
    forma exacta llevando:
      - libres:  un heap con el instante en que cada servidor queda libre
      - en_sist: un heap con los instantes de salida de las solicitudes que
                 siguen en el sistema (para saber cuántas hay al llegar una nueva)

    Devuelve:
      W -> latencia promedio (espera + servicio) de las solicitudes ACEPTADAS
           que llegaron después del calentamiento
      r -> fracción de solicitudes rechazadas después del calentamiento
    """
    libres = [0.0] * c          # todos los servidores libres en t = 0
    heapq.heapify(libres)
    en_sist = []                # salidas pendientes (solicitudes en el sistema)

    t = 0.0
    latencias = []
    llegadas_medidas = 0
    rechazos_medidos = 0

    while True:
        # Siguiente llegada Poisson: tiempo entre llegadas ~ Exp(lam)
        t += rng.exponential(1.0 / lam)
        if t > T_DES:
            break

        # Sacamos del sistema a las solicitudes que ya salieron antes de t
        while en_sist and en_sist[0] <= t:
            heapq.heappop(en_sist)

        medir = t >= WARMUP     # solo contamos estadísticas después del calentamiento
        if medir:
            llegadas_medidas += 1

        # ¿Sistema lleno? (K = capacidad total) -> se rechaza
        if len(en_sist) >= K:
            if medir:
                rechazos_medidos += 1
            continue

        # Aceptada: empieza cuando llegue o cuando se libere el primer servidor
        inicio = max(t, heapq.heappop(libres))
        salida = inicio + rng.exponential(1.0 / mu)
        heapq.heappush(libres, salida)
        heapq.heappush(en_sist, salida)

        if medir:
            latencias.append(salida - t)

    W = np.mean(latencias) if latencias else 1.0 / mu   # sin datos: tiempo de servicio
    r = rechazos_medidos / llegadas_medidas if llegadas_medidas else 0.0
    return W, r


# =============================================================================
# 3. CAPA SD: ECUACIÓN DIFERENCIAL DE U(t)
# =============================================================================
def dU_dt(U, W, r):
    """Tasa de cambio de usuarios activos dados W y r (constantes en el intervalo)."""
    delta = delta0 + kW * max(0.0, W - Wok) + kR * r    # abandono efectivo
    adopcion = (p + q * U / M) * (M - U)                # nuevos usuarios
    return adopcion - delta * U


def rho_pred(U, c):
    """Utilización predicha: carga ofrecida / capacidad total."""
    return a * U / (c * mu)


# =============================================================================
# 4. SIMULACIÓN HÍBRIDA
# =============================================================================
def simular(seed=42):
    rng = np.random.default_rng(seed)
    pasos_por_H = int(round(H / DT))
    n_pasos = int(round(T_FINAL / DT))

    U = float(U0)
    c = c0
    pedido_pendiente = False
    t_llegada_pedido = None
    rho_prev = rho_pred(U, c)

    # Historiales (paso fino del SD)
    t_hist, U_hist, c_hist, rho_hist = [0.0], [U], [c], [rho_prev]
    # Historiales del DES (uno por intercambio)
    tW_hist, W_hist, r_hist = [], [], []
    pedidos, entregas = [], []

    W, r = 0.0, 0.0
    for k in range(n_pasos):
        t = k * DT

        # --- (a) Llegada de servidores pedidos ---
        if pedido_pendiente and t >= t_llegada_pedido - 1e-9:
            c += DELTA_C
            pedido_pendiente = False
            entregas.append(t)
            # Si aun con los nuevos servidores rho_pred > 0.7 -> nuevo pedido inmediato
            if rho_pred(U, c) > RHO_UMBRAL:
                pedido_pendiente = True
                t_llegada_pedido = t + D
                pedidos.append(t)
            rho_prev = rho_pred(U, c)   # la capacidad cambió: actualizamos referencia

        # --- (b) Intercambio SD -> DES cada H días ---
        if k % pasos_por_H == 0:
            lam = a * U                 # llegadas por segundo en la hora pico
            W, r = simular_des(lam, c, rng)
            tW_hist.append(t); W_hist.append(W); r_hist.append(r)

        # --- (c) Integración del SD (Euler) con W y r constantes ---
        U = U + DT * dU_dt(U, W, r)
        U = min(max(U, 0.0), M)         # U no puede salir de [0, M]

        # --- (d) Autoescalado: cruce hacia arriba del umbral ---
        rho = rho_pred(U, c)
        if rho_prev <= RHO_UMBRAL < rho and not pedido_pendiente:
            pedido_pendiente = True
            t_llegada_pedido = t + DT + D
            pedidos.append(t + DT)
        rho_prev = rho

        t_hist.append(t + DT); U_hist.append(U); c_hist.append(c); rho_hist.append(rho)

    return dict(t=np.array(t_hist), U=np.array(U_hist), c=np.array(c_hist),
                rho=np.array(rho_hist), tW=np.array(tW_hist),
                W=np.array(W_hist), r=np.array(r_hist),
                pedidos=pedidos, entregas=entregas)


# =============================================================================
# 5. RESULTADOS Y GRÁFICAS
# =============================================================================
def graficar(res, archivo="resultados_hibrido.png"):
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))

    ax[0, 0].plot(res["t"], res["U"], color="tab:blue")
    ax[0, 0].set_title("Usuarios activos U(t)")
    ax[0, 0].set_xlabel("día"); ax[0, 0].set_ylabel("usuarios")

    ax[0, 1].step(res["t"], res["c"], where="post", color="tab:green")
    for tp in res["pedidos"]:
        ax[0, 1].axvline(tp, color="gray", ls=":", lw=0.8)
    ax[0, 1].set_title("Servidores c(t)  (líneas punteadas = pedidos)")
    ax[0, 1].set_xlabel("día"); ax[0, 1].set_ylabel("servidores")

    ax[1, 0].step(res["tW"], res["W"], where="post", color="tab:orange")
    ax[1, 0].axhline(Wok, color="red", ls="--", label="Wok = 0.75 s")
    ax[1, 0].set_title("Latencia promedio W (DES)")
    ax[1, 0].set_xlabel("día"); ax[1, 0].set_ylabel("segundos"); ax[1, 0].legend()

    ax[1, 1].step(res["tW"], res["r"], where="post", color="tab:red")
    ax[1, 1].set_title("Fracción rechazada r (DES)")
    ax[1, 1].set_xlabel("día"); ax[1, 1].set_ylabel("fracción")

    for a_ in ax.flat:
        a_.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(archivo, dpi=130)
    return fig


if __name__ == "__main__":
    res = simular(seed=42)

    print("=== Resultados (semilla 42) ===")
    print(f"U(100)               = {res['U'][-1]:.1f} usuarios")
    print(f"Servidores finales   = {res['c'][-1]}")
    print(f"Latencia máxima W    = {res['W'].max():.3f} s")
    print(f"Rechazo máximo r     = {res['r'].max():.4f}")
    print(f"Pedidos (día)        = {[round(x, 2) for x in res['pedidos']]}")
    print(f"Entregas (día)       = {[round(x, 2) for x in res['entregas']]}")

    # Robustez: el DES es estocástico, así que repetimos con varias semillas
    print("\n=== 10 réplicas (semillas 0..9) ===")
    filas = []
    for s in range(10):
        rr = simular(seed=s)
        filas.append((rr["U"][-1], rr["c"][-1], rr["W"].max(), rr["r"].max()))
    filas = np.array(filas)
    for nombre, col in zip(["U(100)", "c final", "W max", "r max"], filas.T):
        print(f"{nombre:8s}: media = {col.mean():9.3f}   min = {col.min():9.3f}   max = {col.max():9.3f}")

    graficar(res)
    plt.show()
