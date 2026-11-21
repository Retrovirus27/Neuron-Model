# Estado del proyecto

Última actualización: 2026-11-20

## Hecho
- Exploración inicial del bloque hormonal, en segundos (versión anterior, ya
  no se usa): `wolfram/HormoneModel.wl`, `notebooks/HormoneModel.ipynb`.
  CH-20261120-01.
- Auditoría del modelo publicado (`docs/AUDITORIA_2024.md`). CH-20261120-02.
  Validada por Yordi.
- Memoria del proyecto, plan y skills. CH-20261120-03.
- **Modelo de referencia en horas** con la corrección de POMC:
  `wolfram/HypothalamusModel.wl` + `notebooks/01_modelo_referencia.ipynb`.
  Verificado en Python (V-06). CH-20261120-04.

## Pendiente inmediato
- Yordi: correr `01_modelo_referencia.ipynb` y revisar `validate[]` (V-07, V-08).
- Después: Fase 1 del PLAN (equilibrios, jacobiano y mapa de 24 h de la
  entrada hormonal), dentro del modelo completo de referencia.

## Decisiones abiertas
- ¿Nota de corrección a Sci Rep por las Ec. 7, 9 y 11?
- ¿Herramienta de continuación numérica para la Fase 2?
- Entorno: Wolfram Engine + Jupyter como principal; Wolfram Cloud (gratis)
  solo para exploración interactiva, pegando el código (no permite subir `.nb`).
