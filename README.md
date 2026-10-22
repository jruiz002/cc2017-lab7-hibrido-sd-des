# Laboratorio 7 — Modelo híbrido System Dynamics + simulación de eventos discretos

**Curso:** CC2017 Modelación y Simulación, Universidad del Valle de Guatemala (Sección 10)
**Grupo 5:** q = 0.20, D = 3 días

| Integrante | Carné |
|---|---|
| José Gerardo Ruiz García | 23719 |
| Gerardo André Fernández Cruz | 23763 |
| Juan Diego Solís Martínez | 23720 |

## Modelo

Plataforma de video educativo (Aula Abierta) con dos capas acopladas:

- **Capa SD (días):** usuarios activos `U(t)`,
  `dU/dt = (p + q·U/M)(M − U) − δU`, con `δ = δ0 + κW·max(0, W − Wok) + κR·r`.
  Parámetros: p = 0.002, q = 0.20, M = 20000, U(0) = 100, δ0 = 0.003, κW = 1.0, κR = 0.5, Wok = 0.75 s.
- **Capa DES (hora pico, segundos):** cola M/M/c/K con λ = a·U (a = 0.0005), μ = 2, K = 60.
  Cada intercambio simula 720 s, descarta 120 s de calentamiento y devuelve la latencia media `W`
  y la fracción rechazada `r`.
- **Autoescalado (evento de estado):** `ρ_pred = a·U/(c·μ)`; al cruzar θ = 0.7 hacia arriba sin pedido
  pendiente se piden 2 servidores que llegan D = 3 días después; si al llegar `ρ_pred > θ` se pide otro de inmediato. c(0) = 4.
- **Acoplamiento:** cada H días (base H = 1) el SD entrega λ y c al DES y recibe W y r, que se mantienen
  constantes hasta el siguiente intercambio.

## Estructura del repositorio

```
predicciones/task1_1_predicciones.md   Predicciones escritas antes de simular (no se modifican)
anexo_linea_base_IA/                   Prompt, código y respuesta originales de la IA (Task 1.3, sin modificar)
src/                                   Módulos reutilizables
    des.py          DES M/M/c/K propio y fórmulas de M/M/c (Erlang C, Little)
    sd.py           Capa SD (EDO de Bass con abandono)
    integradores.py Euler y RK4 propios
    eventos.py      Cruce del umbral localizado (bisección y solve_ivp events)
    hibrido.py      Acoplamiento SD-DES con autoescalado (evento de estado y versión ingenua)
notebooks/task1.ipynb, task2.ipynb     Explicación, resultados y gráficas
figuras/                               Gráficas (PNG)
resultados/                            Tablas (CSV)
```

## Cómo reproducir

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
jupyter nbconvert --to notebook --execute --inplace notebooks/task1.ipynb
jupyter nbconvert --to notebook --execute --inplace notebooks/task2.ipynb
```

Todas las simulaciones usan semillas fijas, así que los resultados se repiten exactamente.
