# Conversión de unidades: segundos → horas

Cambio CH-20261120-04. Desde esta versión, el modelo de referencia
(`wolfram/HypothalamusModel.wl`) trabaja con el tiempo en **horas**.
El código original (`reference/NeuronModel.nb`) y `wolfram/HormoneModel.wl`
siguen en segundos (versiones anteriores, solo como referencia; el proyecto
analiza el modelo completo).

| Magnitud | Original (s) | Nuevo (h) | Conversión |
|---|---|---|---|
| Tiempo t | s | h | t_h = t_s / 3600 |
| Ecuaciones de voltaje (V, A, P, L) | C = 1, t en s | `Cm = 1/3600` | dV/dt_h = 3600 · dV/dt_s |
| τ_K | 2 s | 5.556·10⁻⁴ h | τ / 3600 |
| τ_GLU, τ_GABA | 30 s | 8.333·10⁻³ h | τ / 3600 |
| τ_OREXIN | 300 s | 0.0833 h | τ / 3600 |
| τ_increase, τ_decrease (orexina) | 27000, 3312 s | 7.5, 0.92 h | τ / 3600 |
| τ acumulación / eliminación αMSH | 400, 300 s | 0.111, 0.0833 h | τ / 3600 |
| τ acumulación / eliminación AgRP | 300, 300 s | 0.0833, 0.0833 h | τ / 3600 |
| τ activación / inactivación MC4R | 734.4, 2030.4 s | 0.204, 0.564 h | τ / 3600 |
| τ activación / inactivación INSR, GHSR, LEPR | 300 s | 0.0833 h | τ / 3600 |
| k12, k2x, khx, kix | `num[2.30]`, `num[0.71]`, `num[2.354]`, `num[1.03]` | 2.30, 0.71, 2.354, 1.03 h⁻¹ | se usa el valor publicado directamente |
| β, β2 (ghrelina) | `num[1591]`, `num[1238]` | 1591, 1238 h⁻¹ | directo |
| θ, ρ (insulina) | `num[528.63]`, `num[90]` | 528.63, 90 h⁻¹ | directo |
| ε, R0 (leptina) | `num[6231.39]`, `num[0.11]` | 6231.39, 0.11 h⁻¹ | directo |
| Ingesta (tasa durante la comida) | `num[2.1]`, `num[1.1]`, `num[0.8]` | 2.1, 1.1, 0.8 L/h | directo |
| Duración de desayuno / almuerzo / cena | 2340 / 1676 / 1400 s | 0.650 / 0.466 / 0.389 h | / 3600 |
| ω_c (circadiano) | 0.000073 rad/s | 0.2628 rad/h | × 3600 |
| Conductancias, potenciales, pendientes, EC50, concentraciones, σ, γ, δ, s, η, l, FM, BV | – | – | sin cambio |

## Cómo se validó
1. Port en Python en horas (`analysis/reference_model/model_h.py`), escrito con
   los parámetros convertidos (no reescalando la función en segundos).
2. Comparado con la versión en segundos corregida (`analysis/audit_2024/`):
   - Hormonas: diferencia relativa máxima de 5·10⁻⁷ (TRE) y 1·10⁻⁵ (CD).
   - Spikes por bin de 30 min: diferencia máxima de ±4; totales dentro del
     1 %. Es del mismo orden que la diferencia puramente numérica entre
     dt = 0.01 s y 0.005 s.
3. La referencia para Wolfram son los CSV en `analysis/reference_model/ref/`
   (modelo corregido en horas, 24 h, dietas TRE, CD, SM e IF).
