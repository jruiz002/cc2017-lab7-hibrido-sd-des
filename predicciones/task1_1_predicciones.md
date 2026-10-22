# Task 1.1 — Predicciones (grupo 5: q = 0.20, D = 3 días)

Con g = 5: (5−1) mod 4 = 0, así que q = 0.20 y D = 3 días.

a. Retroalimentación completa vs. cascada

U(100): menor con el lazo cerrado, aunque la diferencia debería ser pequeña (estimamos entre 2% y 8% menos). Con el lazo cerrado, la latencia alta antes de que lleguen servidores provoca abandono y frena el crecimiento. Cuando llegan los servidores, la latencia baja y la plataforma se recupera casi por completo.
Latencia máxima: menor con el lazo cerrado. En la cascada los usuarios siguen creciendo aunque la cola esté saturada. En el lazo, el abandono reduce la carga, así que funciona como un freno natural.
Fracción máxima rechazada: cercana a 0 en ambos casos, un poco menor en el lazo. Con K = 60 y ρ < 1 la cola casi nunca se llena; el riesgo solo existe si U pasa de 16,000 con c = 4.

b. Primer pedido y servidores al día 100

El umbral con c = 4 es aU/8 = 0.7, es decir U ≈ 11,200. Con crecimiento tipo Bass a una tasa efectiva de ~0.2 por día desde U(0) = 100, estimamos el primer pedido alrededor del día 26–30, con llegada 3 días después (c = 6).
Con c = 6 el siguiente umbral es U ≈ 16,800. Como sin congestión el techo es ≈ 20,000·(1 − δ0/q) ≈ 19,700, sí debería cruzarse, alrededor del día 40–50, y c pasaría a 8.
Con c = 8 el umbral sería U = 22,400 > M, así que no habría más pedidos. Predicción: 8 servidores al día 100 (6 si el abandono por latencia retrasa demasiado el crecimiento).

c. H = 1 vs. H = 7

La métrica que más cambia es la latencia máxima, que aumenta. Con H = 7 el SD usa una W “vieja” durante una semana: no se entera a tiempo de la congestión ni de que ya llegaron los servidores, y el abandono responde tarde.
U(100) debería cambiar poco, porque al día 100 el sistema ya está cerca del equilibrio.
Esperamos un cambio irregular con H, porque depende de cómo se alinean los intercambios con la llegada de servidores (D = 3 es menor que H = 7) y del ruido del DES en cada ventana.
