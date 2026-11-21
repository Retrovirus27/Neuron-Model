# Plan de trabajo: análisis de sistemas dinámicos del modelo hipotalámico

Objeto de estudio: el **modelo completo** (hormonas, receptores y red de
POMC, AgRP, PVH y LHA). El bloque hormonal se analiza como la entrada que
fuerza a la red, no como un modelo aparte.

Base: Tamayo-Molina et al. (2024) *Sci Rep* 14:28897. Fecha del plan: 2026-11-20.

## Diagnóstico de partida (lo que el nuevo trabajo debe resolver)
1. **El "sistema oscilatorio doble" no es endógeno.** No hay ningún oscilador
   autónomo: LHA oscila porque la fuerza `ICircadian`, y POMC/AgRP porque los
   fuerzan las comidas. Es un sistema forzado por dos entradas periódicas. Para
   hablar de un *food clock* endógeno hace falta un ciclo límite que persista
   o que se arrastre al patrón de comidas.
2. **Validación neuronal cualitativa.** Las hormonas se comparan con datos;
   la actividad neuronal se interpreta a ojo. Faltan métricas (fase, amplitud,
   desfase respecto al ritmo circadiano).
3. **No hay análisis de sensibilidad ni de robustez** (parámetros de fuentes
   y especies distintas).
4. **Discrepancias paper ↔ código** (ver `AUDITORIA_2024.md`): dos son errores
   de transcripción del texto (leptina, ghrelina) y una es un error real del
   código (I_LEPR en POMC con el voltaje de AgRP).
5. **El subsistema hormonal es una cascada sin retroalimentación:** para cada
   ingesta constante tiene un único equilibrio estable (no hay multiestabilidad
   posible). Con comidas periódicas converge a una única órbita de 24 h.

## Fase 0. Base limpia
- [x] Exploración inicial del bloque hormonal (`wolfram/HormoneModel.wl`;
      versión anterior, reemplazada por el modelo de referencia).
- [x] Auditoría del paper vs. código vs. datos (`AUDITORIA_2024.md`).
- [x] Versión de referencia del modelo completo en Wolfram con I_LEPR
      corregido, en un solo `modelEqs[diet]` (CH-20261120-04; falta V-07).
- [x] Horarios periódicos con `Mod[t, 24 h]` (comidas y β de la ghrelina) para
      simular cualquier número de días; IFI sólo cubría el día 1.
- [x] Unidades: horas (`docs/UNIDADES.md`). Adimensionalizar queda para la
      Fase 2 si hace falta.
- [ ] Comidas suaves (pulsos sigmoides) además de los escalones, para poder
      hacer continuación numérica.
- [x] Revisar nombres que chocan con Wolfram (`I`, `h[h_,…]`).

## Fase 1. La entrada hormonal del modelo completo (rápido de hacer, base del resto)
Se caracteriza con precisión la señal que fuerza a la red, para usarla
después en las Fases 2 y 3.
- [ ] Equilibrios analíticos en función de la ingesta constante F (`Solve`).
- [ ] Jacobiano simbólico, valores propios → tiempos de relajación de cada
      hormona.
- [ ] Respuesta periódica a las comidas: mapa estroboscópico de 24 h,
      convergencia a la órbita periódica, amplitudes y fases según la dieta.
- [ ] Sensibilidad local (derivadas) y global (Sobol) de amplitud/fase.

## Fase 2. Análisis rápido–lento (aporte central)
Escalas: rápida (spikes), intermedia (sinapsis, orexina τ ≈ 7.5 h),
lenta (hormonas, horas).
- [ ] Congelar a_INSR, a_GHSR, a_LEPR e I_circadian como **parámetros**.
- [ ] Diagrama de bifurcación de cada población aislada y de la red:
      reposo ↔ disparo tónico ↔ bursting (silla-nodo, Hopf, homoclínica).
- [ ] Curvas de frecuencia de disparo vs. activación de receptor (f–a).
- [ ] Superponer la trayectoria lenta diaria de las hormonas sobre el diagrama:
      explicar mecánicamente cuándo se encienden/apagan POMC y AgRP
      (p. ej. "comidas pequeñas no alcanzan el umbral" = la insulina no cruza
      el punto de bifurcación).
- [ ] Reducción a modelo de tasa de disparo (*firing-rate/mean-field*) para
      que la red completa sea tratable.

## Fase 3. Cuantificar alineación / desalineación
- [ ] Definir la fase de cada población (centro de masa de su actividad diaria)
      y el desfase respecto a `ICircadian`.
- [ ] Barrer horario y frecuencia de comidas → mapas de alineación.
- [ ] Escenarios de jet lag y comer por turnos (desplazar la fase de comidas).

## Fase 4. Cerrar lazos (aquí aparecen varios estados estables)
- [ ] **Ingesta controlada por las neuronas** (AgRP/PVH disparan la comida):
      el *food clock* se vuelve endógeno (ciclo límite propio); estudiar
      entrainment a la luz, lenguas de Arnold, jet lag.
- [ ] **Lazo lento leptina ↔ masa grasa** (FM dinámico, balance energético
      tipo Tam 2009): buscar biestabilidad e histéresis (peso normal vs.
      obesidad).
- [ ] Revisar si AgRP ⊣ POMC con LHA forma un *flip-flop* (coexistencia de
      dos estados en algún rango de parámetros).

## Validación
Cuando el modelo se corrige a propósito (p. ej. la corriente de POMC), la
Fig. 2A original deja de ser la referencia exacta. La referencia es el port
en Python del modelo vigente (`analysis/reference_model/`), verificado contra
el código publicado cambiando una sola cosa a la vez.

## Herramientas
- Wolfram (principal): `Solve`, `D`, `Eigenvalues`, `ParametricNDSolve`,
  `StreamPlot`, barridos. En Jupyter no funciona `Manipulate`
  (usar Wolfram Cloud o Mathematica para exploración interactiva).
- Continuación de ramas (Hopf, silla-nodo, ciclos): XPPAUT/AUTO, MatCont o
  BifurcationKit.jl (complemento, decidir al llegar a la Fase 2).
- Verificación cruzada: port en Python (`analysis/audit_2024/model.py`).

## Orden recomendado
Fase 0 → Fase 1 → Fase 2 (AgRP y POMC vs. activación de receptores).
Con eso ya hay un paper sólido; la Fase 4 puede ser el gancho principal o un
tercer trabajo.
