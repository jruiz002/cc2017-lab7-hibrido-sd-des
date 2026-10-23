Hola, estamos haciendo un laboratorio del curso de Modelación y Simulación y necesitamos que nos ayudes a programar en Python un modelo híbrido que combina System Dynamics con simulación de eventos discretos. Te pasamos el modelo completo:
Capa SD (tiempo en días): U(t) son los usuarios activos.
dU/dt = (p + q*U/M)*(M - U) - delta*U
delta = delta0 + kW*max(0, W - Wok) + kR*r
Parámetros: p = 0.002, q = 0.20, M = 20000, U(0) = 100, delta0 = 0.003, kW = 1.0, kR = 0.5, Wok = 0.75 segundos.
Capa DES (hora pico): cola M/M/c/K. Llegadas Poisson con lambda = a*U, a = 0.0005 solicitudes por segundo por usuario. Servicio exponencial con media 0.5 s (mu = 2). Capacidad máxima K = 60, si está lleno se rechaza la solicitud. Cada vez que el SD necesita datos, el DES simula 720 segundos, descarta los primeros 120 como calentamiento y devuelve la latencia promedio W de las solicitudes aceptadas y la fracción rechazada r.
Autoescalado: rho_pred = a*U/(c*mu). Cuando rho_pred pasa hacia arriba de 0.7 y no hay pedido pendiente, se piden 2 servidores que llegan D = 3 días después. Si cuando llegan rho_pred sigue arriba de 0.7, se pide otro pedido de inmediato. Al inicio c = 4.
Acoplamiento: cada H días el SD le pasa al DES lambda y c, y el DES devuelve W y r, que se mantienen constantes hasta el siguiente intercambio. H = 1 día.
Necesitamos que el código simule 100 días y nos muestre U(100), el número final de servidores, la latencia máxima y la fracción máxima de rechazos, con gráficas de U(t), c(t), W y r. Por favor coméntalo para entenderlo y dinos si los resultados tienen sentido.
