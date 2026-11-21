---
name: analisis-dinamico
description: Flujo de trabajo para el análisis de sistemas dinámicos del modelo hipotalámico (equilibrios, jacobiano y estabilidad, análisis rápido–lento, bifurcaciones, órbitas periódicas forzadas, entrainment, multiestabilidad). Úsala al planear o ejecutar cualquier fase de docs/PLAN.md.
---

# Análisis de sistemas dinámicos del modelo

Leer antes `docs/PLAN.md` (fases) y `docs/ESTADO.md` (qué está hecho).
El código va en Wolfram (ver skill `wolfram-jupyter`).

## Hechos estructurales que no hay que redescubrir
- El subsistema hormonal (GL1, GL2, ghrelina, insulina, leptina) es una
  **cascada sin retroalimentación**. Con ingesta constante F tiene un único
  equilibrio (GL2* = F/k2x), y el jacobiano es casi triangular con valores
  propios negativos (−k_hx, −k_ix, −ε/BV y los del bloque de glucosa). No puede
  haber multiestabilidad hormonal. Equilibrio en ayuno: ghrelina = β/k_hx ≈ 676,
  insulina = ρ/k_ix ≈ 87, leptina ≈ 22 (unidades del modelo).
- El sistema completo es **forzado** (comidas + `ICircadian`), no autónomo:
  sus atractores son órbitas periódicas de 24 h, no equilibrios. Para hablar
  de "estados estables" usar el mapa estroboscópico (muestrear cada 24 h).
- Escalas de tiempo: spikes (s), sinapsis (30–400 s), orexina
  (τ_inc = 27000 s), receptores hormonales (300 s), hormonas (horas). Esto
  justifica el análisis rápido–lento.

## Recetas en Wolfram
- **Equilibrios:** `Solve[rhs == 0 /. food -> F, vars]` o `NSolve`;
  comprobar unicidad.
- **Jacobiano:** `J = D[rhs, {vars}]`; estabilidad con `Eigenvalues[J /. eq]`;
  tiempos de relajación = −1/Re(λ).
- **Barridos:** `ParametricNDSolveValue` + `Table`; ver los atractores
  descartando un transitorio (≥ 3–5 días en el modelo completo).
- **Rápido–lento:** congelar a_INSR, a_GHSR, a_LEPR e I_circadian como
  parámetros; para cada valor, simular la población (o la red) y medir:
  reposo, tónico o bursting; frecuencia; voltaje mín/máx. Graficar el
  diagrama y superponer la trayectoria lenta diaria.
- **Fase / alineación:** fase de una población = centro de masa circular de
  su actividad diaria (`Arg[Total[r Exp[2 Pi I t/24]]]`); desfase respecto
  al máximo de `ICircadian`.
- **Continuación de ramas** (Hopf, silla-nodo, ciclos): Wolfram no la trae.
  Opciones: XPPAUT/AUTO, MatCont o BifurcationKit.jl. Requiere forzamiento
  suave (comidas sigmoides) o el sistema autónomo congelado.

## Para multiestabilidad hay que cerrar lazos
- Ingesta controlada por la actividad de AgRP/PVH (oscilador endógeno).
- Masa grasa dinámica regulada por la leptina (biestabilidad del peso).
- Inhibición mutua tipo *flip-flop* en la red.
Cada lazo nuevo es un cambio de modelo: registrar, justificar biológicamente
y validar que, con el lazo desactivado, se reproduce la Fig. 2A.

## Verificación
- Cada resultado numérico se contrasta con una segunda vía (analítica vs.
  numérica, o Wolfram vs. el port en Python).
- Revisar sensibilidad al paso (`MaxStepSize`) y a la duración del transitorio.
- Nada se declara validado sin confirmación de la persona usuaria.
