# Auditoría del modelo publicado (Sci Rep 2024)

Fecha: 2026-11-20 (CH-20261120-02). Método: port del notebook
`reference/NeuronModel.nb` a Python (`analysis/audit_2024/model.py`, RK4,
dt = 0.01 s, spikes = cruce de −20 mV, bins de 30 min). La versión
"publicada" reproduce la Fig. 2A del paper, así que el port es fiel.

Para regenerar: dentro de `analysis/audit_2024/` correr
`python3 runall.py TRE 24 publicado,LEPR_POMC_corregido,leptina_como_texto,ghrelina_como_texto`
(igual con `CD 24` y con `TRE 72` / `CD 72`) y luego `python3 figs.py`.
Requiere `numpy`, `numba`, `matplotlib`.

## Discrepancias encontradas

| # | Dónde | Texto del paper | Código (figuras) | Tipo |
|---|---|---|---|---|
| 1 | Leptina, Ec. 9 | `R = R0 η/(1+(IN/EC50)^(−l))`: la insulina estimula | `^(+l)`: la insulina inhibe | Error de transcripción del texto |
| 2 | Ghrelina, Ec. 7 | β1 en 0–6 h, β2 en 6–24 h | β2 (más bajo) en 3–6 h, β el resto | Error de descripción del texto |
| 3 | POMC, Ec. 1 | I_LEPR sobre POMC | `ILEPR[A[t],…]`: usa el voltaje de AgRP | **Error del código** |
| 4 | I_circadian, Ec. 11 | fase + π/3 | fase − π/3 | Error de transcripción del texto (las figuras usan − π/3) |
| 5 | Dieta IF (Fig. 3) | – | `IFI` solo define comidas el día 1 | Limitación del código (las figuras de IF son de 24 h, no se afectan) |

## Impacto

### 1 y 2: el texto, no los resultados
Las figuras salen del código, y el código es la forma ajustada a los datos
humanos (Cummings 2001). Si se implementa el modelo tal como está escrito en el
texto, con los mismos parámetros:

| Versión | RMSE ghrelina (pM) | RMSE insulina (pM) | RMSE leptina (pM) |
|---|---|---|---|
| Código publicado | 22.9 | 88.5 | 0.78 |
| Leptina como Ec. 9 | 22.9 | 88.5 | **6.97** (cae a ~2–4 pM) |
| Ghrelina como texto | **41.0** | 88.5 | 0.78 |

Las conclusiones del paper no cambian. El problema es de **reproducibilidad**:
alguien que reimplemente desde el texto no obtiene las figuras. El código está
en Mendeley, y eso mitiga el problema.

### 3: efecto real pero acotado
Comparación publicado → corregido (total de spikes y correlación entre las
series de 30 min):

| Dieta | AgRP | POMC | PVH | LHA |
|---|---|---|---|---|
| 3 comidas, 24 h | −5 %, r = 0.97 | −10 %, r = 0.89 | −9 %, r = 0.88 | −1 %, r = 0.87 |
| Dieta común, 24 h | −10 %, r = 0.96 | −12 %, r = 0.88 | −12 %, r = 0.88 | −1 %, r = 0.80 |
| TRE, 72 h | −6 % | −9 % | −7 % | −2 % |
| Dieta común, 72 h | −17 % | −12 % | −13 % | −1 % |

- **Se mantienen** las conclusiones principales: actividad nocturna de
  POMC/PVH/AgRP, LHA diurna (~17 h), POMC/PVH siguen a la insulina, TRE
  alineado durante 3 días y dieta común con pérdida de la actividad nocturna
  de POMC/PVH.
- **Cambian detalles que el texto interpreta:**
  - La alternancia POMC/PVH vs. LHA entre 2 y 4 h (la "fragmentación del
    sueño") **se mantiene**; con la corrección los picos quedan más definidos
    (revisado por Yordi, 2026-11-20).
  - En la dieta común, la pequeña reactivación de POMC/PVH tras la merienda de
    las 16 h y los repuntes de AgRP en la tarde/noche desaparecen.
  - POMC nocturno baja de ~41 a ~33 spikes/30 min.

Figuras: `analysis/audit_2024/figures/fig1_hormonas.png`,
`fig2_neuronas_24h.png`, `fig3_72h.png`.

## Decisión
- Mantener las formas ajustadas del código para la leptina y la ghrelina.
- Usar I_LEPR corregido en todo el trabajo nuevo.
- Para el paper nuevo: mencionar la corrección de forma transparente (en
  Métodos o en material suplementario). Valorar con la persona usuaria si se
  envía una nota de corrección (*correction*) a Sci Rep por el texto de las
  Ec. 7 y 9.
