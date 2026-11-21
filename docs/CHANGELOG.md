# CHANGELOG

## CH-20261120-05: preparación del repositorio
- `.gitignore` (checkpoints de Jupyter, caches de Python, `.pkl`, `.DS_Store`, `backup/`).
- Rama de trabajo propuesta: `dynamical-systems` en
  github.com/Retrovirus27/Neuron-Model. `main` queda intacta con el código
  del paper.

## CH-20261120-04: modelo de referencia en horas
- `wolfram/HypothalamusModel.wl`: modelo completo con tiempo en horas
  (tabla en `docs/UNIDADES.md` y `unitTable`), corriente de leptina en POMC
  corregida, una sola `modelEqs[diet]`, horarios periódicos con `Mod[t, 24]`
  (comidas y β de la ghrelina), nombres seguros, detección de spikes con
  `WhenEvent`, `simulateHormones` (solo hormonas), `plotSim` y `validate`.
- `notebooks/01_modelo_referencia.ipynb`: carga, simula las 4 dietas y valida.
- Referencia Python en horas: `analysis/reference_model/` (`model_h.py`,
  `validate_h.py`, `refcsv.py`, `ref/*.csv`). `validate_h.py` necesita los
  `.pkl` de `analysis/audit_2024` (correr antes `runall.py`).
- Consecuencia: IF ahora tiene comidas todos los días (en el original solo
  el día 1).
- `wolfram/HormoneModel.wl` y `notebooks/HormoneModel.ipynb` quedan como
  versión anterior (en segundos); no se borran.

## CH-20261120-03: memoria del proyecto, plan y skills
- `AGENTS.md`, `CLAUDE.md`, `docs/PLAN.md`, `docs/ESTADO.md`,
  `docs/VALIDACION.md`, `docs/Promps.md`.
- Skills en `skills/` (copiar a `.claude/skills/` para Claude Code): `wolfram-jupyter`, `auditoria-modelo`,
  `analisis-dinamico`.
- Copias de referencia en `reference/` (paper y `NeuronModel.nb`).

## CH-20261120-02: auditoría del modelo publicado
- Port del modelo completo (25 ODE) a Python + numba:
  `analysis/audit_2024/model.py`, `runall.py`, `figs.py`.
- Reproduce la Fig. 2A del paper.
- Compara 4 variantes: publicado, I_LEPR corregido, leptina como texto,
  ghrelina como texto. Resultados en `docs/AUDITORIA_2024.md`.

## CH-20261120-01: modelo hormonal desacoplado
- `wolfram/HormoneModel.wl` y `notebooks/HormoneModel.ipynb`: GL1, GL2,
  ghrelina, insulina, leptina; dietas TRE, CD, IF; gráficas 24/72 h con datos
  de Cummings 2001; anexo con activación de receptores.
- Cambios vs. original: nombres ASCII para griegas, `i[gl2_, I_]` → `dINS`,
  `MaxStepSize -> 60` s en `NDSolve` (sin neuronas el solver podría saltarse
  las comidas).
