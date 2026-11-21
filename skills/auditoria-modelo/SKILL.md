---
name: auditoria-modelo
description: Auditar el modelo hipotalámico comparando paper, código y datos ajustados, y medir el impacto de cada discrepancia con simulaciones antes y después. Úsala cuando aparezca una diferencia entre ecuaciones del texto y el código, o antes de cambiar una ecuación del modelo.
---

# Auditoría paper ↔ código ↔ datos

Referencias: `reference/TamayoMolina2024_SciRep.pdf`,
`reference/NeuronModel.nb` (intocables) y `docs/AUDITORIA_2024.md`
(auditoría previa y decisiones).

## Principio
**Los datos ajustados mandan.** Los parámetros hormonales se ajustaron a datos
humanos (Cummings 2001) con la forma escrita en el código. Una discrepancia
con el texto no significa que el código esté mal: primero se mide cuál de las
dos formas reproduce los datos y las figuras publicadas.

## Procedimiento
1. **Extraer el código.** Para leer el `.nb` sin Wolfram, parsear las celdas
   `Input` (ver el parser usado en la sesión de 2026-11-20: convierte
   `RowBox`/`FractionBox`/`SubscriptBox` a texto). Ojo con el orden de
   evaluación: hay celdas que usan parámetros definidos más abajo con `:=`.
2. **Tabular cada ecuación**: forma en el texto, forma en el código y
   parámetros. Clasificar cada diferencia:
   - *Transcripción del texto*: el código reproduce las figuras y los datos;
     el texto no.
   - *Error del código*: contradice la descripción biológica o el circuito
     (p. ej. una corriente calculada con el voltaje de otra población).
   - *Elección de modelado*: decisión defendible que conviene justificar.
3. **Simular las variantes** con `analysis/audit_2024/model.py`. Las variantes
   se activan con `flags` en `rhs`; para una variante nueva, agregar un flag sin
   cambiar la versión publicada. Usar RK4 a dt = 0.01 s y verificar la
   convergencia con dt = 0.005 s.
4. **Medir:**
   - Hormonas: RMSE frente a los datos de Cummings (ghrelina ÷3.3709 y
     leptina ÷2.0043 para pasar a pM).
   - Neuronas: spikes por bin de 30 min (cruce de −20 mV), total de spikes,
     correlación entre versiones, episodios de actividad.
   - Revisar cada afirmación del texto que dependa del cambio.
5. **Primero reproducir**: la variante "publicado" debe reproducir la Fig. 2A
   antes de interpretar diferencias.
6. **Registrar**: actualizar `docs/AUDITORIA_2024.md` (o un documento nuevo
   de auditoría), `docs/CHANGELOG.md` y `docs/VALIDACION.md`, con figuras en
   `analysis/<auditoria>/figures/`.

## Comunicación
Ser directo sobre qué cambia y qué no. Distinguir error de reproducibilidad
(texto) de error de resultado (código). No dramatizar ni minimizar: dar
números.
