# Validación

| Caso | Cambio | Cómo se valida | Estado |
|---|---|---|---|
| V-01 | CH-20261120-01 | Rangos de las 5 ODE en Python: insulina ~87–460 pM, ghrelina ~110–200 pM, leptina ~8–11 pM (día 2, TRE) | Verificado por Claude (Python) |
| V-02 | CH-20261120-01 | Correr `notebooks/HormoneModel.ipynb` en Jupyter/Wolfram; las curvas deben coincidir con la Fig. 2A | Reemplazado por V-07 (versión anterior) |
| V-03 | CH-20261120-02 | La variante "publicado" del port reproduce la Fig. 2A (AgRP ~50 hasta las 8 h, POMC/PVH tras cada comida, LHA ~17 h activa) | Verificado por Claude (comparación visual con el paper) |
| V-04 | CH-20261120-02 | Convergencia numérica: dt 0.01 vs 0.005 s → spikes iguales (±2) | Verificado por Claude |
| V-05 | CH-20261120-02 | Revisión de la auditoría y sus conclusiones | Validado por Yordi (2026-11-20): la fragmentación y las oscilaciones se mantienen; los picos quedan más definidos |
| V-06 | CH-20261120-04 | Port en horas vs. versión en segundos corregida (Python): hormonas ≤ 1e-5 relativo; spikes ±4 por bin, totales dentro del 1 % | Verificado por Claude |
| V-07 | CH-20261120-04 | `notebooks/01_modelo_referencia.ipynb`: `validate[]` con resultado OK en TRE, CD, SM, IF (hormonas < 1e-3; spikes ±3 % del total y ≤ 5 por bin) | **Pendiente (Yordi)** |
| V-08 | CH-20261120-04 | Figuras de `plotSim` coherentes con la Fig. 2A, salvo los cambios esperados por la corrección de POMC (picos más definidos) | **Pendiente (Yordi)** |
