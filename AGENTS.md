# Neuron_model: memoria operativa para sesiones futuras

Antes de modificar nada, lee `docs/ESTADO.md` y las últimas entradas de
`docs/CHANGELOG.md`. Para cambios de alcance o dirección, lee `docs/PLAN.md`.
Para cualquier dato del modelo publicado, la referencia es `reference/` y
`docs/AUDITORIA_2024.md`.

## Qué es
Continuación del modelo del sistema hipotalámico publicado en
Tamayo-Molina et al., *Sci Rep* 14:28897 (2024)
(`reference/TamayoMolina2024_SciRep.pdf`, código original en
`reference/NeuronModel.nb`). El modelo une hormonas en sangre (glucosa
GI1/GI2, ghrelina, insulina, leptina), su detección por receptores
(INSR, GHSR, LEPR) y una red de 4 poblaciones neuronales tipo conductancia
(POMC, AgRP, PVH, LHA) forzada por un ritmo circadiano.

**Meta del nuevo trabajo:** análisis de sistemas dinámicos del **modelo
completo** (hormonas + receptores + red neuronal), no solo del bloque
hormonal: equilibrios, estabilidad, análisis rápido–lento, bifurcaciones, entrainment
y, si se cierran lazos, multiestabilidad). Ver `docs/PLAN.md`.

Trabajo entre la persona usuaria (Yordi) y Claude. No hay más coautores en
esta fase.

## Estructura
- `reference/`: paper y notebook original. **Intocable**: no se modifica,
  mueve ni borra.
- `wolfram/`: código Wolfram Language (`.wl`) que se carga con `Get[]`.
- `notebooks/`: notebooks `.ipynb` para JupyterLab con kernel Wolfram.
- `analysis/`: análisis auxiliares (p. ej. `audit_2024/`, port en Python).
- `docs/`: memoria del proyecto (PLAN, ESTADO, CHANGELOG, VALIDACION,
  AUDITORIA_2024, Promps).
- `skills/`: skills del proyecto (ver abajo). Para que Claude Code las cargue solas, copiarlas a `.claude/skills/` (`cp -r skills .claude/`).

## Límites del proyecto
- **Lenguaje principal: Wolfram Language** (experiencia de la persona
  usuaria). Se ejecuta en JupyterLab con el kernel Wolfram (Wolfram Engine 15,
  entorno conda `wolfram`). Python solo para verificaciones cruzadas o
  herramientas que Wolfram no tenga (p. ej. continuación numérica).
- Repositorio: github.com/Retrovirus27/Neuron-Model. `main` tiene el código
  del paper 2024 y no se modifica. Rama de trabajo: `dynamical-systems`.
- **Git: nunca hacer commit ni push.** Al cerrar un bloque, proponer los
  mensajes de commit (en inglés: `model: ...`, `analysis: ...`, `docs: ...`).
  La persona usuaria los ejecuta.
- **Nunca borrar archivos** sin indicación explícita. Solo crear o modificar.
- No editar un archivo mientras la persona usuaria lo está ejecutando; si solo
  pide una respuesta, contestar en texto sin editar.
- En el entorno de Claude **no hay Wolfram**: el código Wolfram solo se revisa;
  la validación numérica la hace la persona usuaria en su Mac. Claude puede
  verificar cualquier resultado con el port en Python
  (`analysis/audit_2024/model.py`).
- Entorno: conda (no mamba).

## No romper (decisiones tomadas)
- **Los parámetros hormonales están ajustados a datos humanos**
  (Cummings et al. 2001). La forma que manda es la del **código**, no la del
  texto del paper:
  - Leptina: `R0 η / (1 + (IN/EC50)^(+l))` (la insulina inhibe la producción
    de leptina en la escala de horas). La Ec. 9 del paper escribe `^(-l)`:
    es un error de transcripción. Con `^(-l)` y los mismos parámetros el
    ajuste se rompe (RMSE leptina 0.78 → 6.97 pM).
  - Ghrelina: `β2` (más bajo) entre las 3 y las 6 h de cada día y `β` el resto.
    El texto dice β1 en 0–6 h y β2 en 6–24 h: error de descripción
    (RMSE ghrelina 22.9 → 41.0 pM).
  - No "corregir" estas formas para que coincidan con el texto.
- **Error real del código publicado:** en la ecuación de POMC, la corriente de
  leptina usa el voltaje de AgRP: `ILEPR[A[t],ALEPR[t],-20]`. Lo correcto es
  `ILEPR[P[t],ALEPR[t],-20]`. En el trabajo nuevo se usa la versión corregida
  (impacto en `docs/AUDITORIA_2024.md`).
- Unidades del modelo original: tiempo en segundos; `num[x] = x/3600` pasa
  tasas por hora a por segundo; `gen[h]` pasa horas a segundos. Cualquier
  cambio de unidades (p. ej. trabajar en horas) se registra en CHANGELOG
  y se valida reproduciendo la Fig. 2A.
- Objeto de estudio: el modelo completo (`wolfram/HypothalamusModel.wl`).
  El bloque hormonal no depende de las neuronas (acoplamiento en un solo
  sentido), así que funciona como **entrada** del sistema y se puede resolver
  rápido con `simulateHormones`, pero el análisis es del sistema completo.
  `wolfram/HormoneModel.wl` y `notebooks/HormoneModel.ipynb` son una versión
  anterior (en segundos) y ya no se usan.
- Reproducción de referencia: el port en Python reproduce la Fig. 2A del paper
  (dt = 0.01 s, RK4; spikes = cruce de −20 mV hacia arriba, bins de 30 min).
  Cualquier versión nueva del modelo debe reproducir la Fig. 2A antes de
  usarse para análisis.

## Registro obligatorio
Cada cambio se registra en `docs/CHANGELOG.md` con ID `CH-YYYYMMDD-NN` y sus
casos en `docs/VALIDACION.md`. `docs/ESTADO.md` se actualiza al terminar.
También se incluyen los prompts utilizados en `docs/Promps.md`.
Excepción: cuando la persona usuaria diga **"No registrar este prompt"** (o lo
pida de forma equivalente), ese prompt y lo que se haga a partir de él no se
registran en `docs/CHANGELOG.md` ni en `docs/Promps.md`. Los demás documentos
(ESTADO, PLAN, AGENTS) sí se corrigen si el cambio lo requiere.
Nada se declara validado sin confirmación de la persona usuaria (que corre el
código Wolfram y avisa si algo falla: entonces se corrige o se revierte).

## Skills del proyecto
- `wolfram-jupyter`: cómo escribir y entregar código Wolfram para el kernel
  de Jupyter.
- `auditoria-modelo`: revisar consistencia paper ↔ código ↔ datos ajustados
  y medir el impacto de cada discrepancia.
- `analisis-dinamico`: flujo de trabajo para equilibrios, estabilidad,
  rápido–lento, bifurcaciones y entrainment.

## Co-autoría
En los commits solo incluir `Co-Authored-By: Claude Opus 5.5`. No incluir el
email ni la sesión.
