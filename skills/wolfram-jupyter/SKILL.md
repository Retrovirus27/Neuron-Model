---
name: wolfram-jupyter
description: Escribir, revisar o entregar código Wolfram Language para este proyecto (archivos .wl y notebooks .ipynb con kernel Wolfram en JupyterLab). Úsala antes de crear o editar cualquier código del modelo.
---

# Código Wolfram para JupyterLab (kernel Wolfram)

Entorno de la persona usuaria: macOS, Wolfram Engine 15 (gratis), JupyterLab
4.x en el entorno conda `wolfram`, kernel registrado con
WolframLanguageForJupyter. En el entorno de Claude **no hay Wolfram**: el
código no se puede ejecutar aquí. Hay que revisarlo con cuidado y, cuando haya
números, verificarlos con el port en Python (`analysis/audit_2024/model.py`)
antes de entregar.

## Dónde va cada cosa
- Definiciones reutilizables (parámetros, ecuaciones, funciones) → `wolfram/*.wl`.
  Se cargan con `Get["<ruta>/archivo.wl"]`.
- Exploración, figuras y pasos del análisis → `notebooks/*.ipynb`.
  Una ecuación o idea por celda, con su fórmula en LaTeX en una celda markdown
  arriba.
- Generar el `.ipynb` por script (JSON nbformat 4), con
  `kernelspec.name = "wolframlanguage15"`. Si el nombre no coincide, Jupyter
  pide elegir kernel: avisar a la persona usuaria.

## Reglas de escritura
- Primera celda: `ClearAll["Global`*"]`.
- Nombres en ASCII (`sigma`, `beta`…): se teclean fácil en Jupyter.
  Documentar el mapeo al símbolo original.
- No usar como variable ni como patrón símbolos del sistema o muy cortos y
  peligrosos: `I`, `E`, `C`, `D`, `N`, `K`, `O`. No reutilizar el nombre de
  la función como argumento (`h[h_, …]`).
- Mayúsculas importan en estilos y opciones (`"Input"`, no `"input"`).
- Respetar las unidades del modelo y documentarlas en la celda: en el original
  el tiempo va en segundos (`gen[h]` = h·3600, `num[x]` = x/3600).
- `NDSolve` con entradas discontinuas (comidas tipo escalón): fijar
  `MaxStepSize` menor que la comida más corta (≈ 60 s), o usar
  `WhenEvent`/pulsos suaves.
- Sistemas grandes: usar `NDSolveValue` y graficar con `Plot` sobre la
  interpolación. Evitar tablas a 0.1 s salvo que se necesiten (conteo de
  spikes).
- Opciones en una lista dentro de funciones `HoldAll` (`Plot`):
  `Evaluate[Sequence @@ opts]`.
- `Manipulate`, `Dynamic` y los sliders **no funcionan** en Jupyter. Para
  barridos usar `Table`/`ParametricNDSolve` + `ListPlot`/`GraphicsGrid`.
  Para exploración interactiva, sugerir Wolfram Cloud.
- La salida larga se trunca en Jupyter: para inspeccionar, usar
  `Scan[Print, …]` o exportar a archivo.
- Figuras: `Export["figures/<nombre>.png", fig, ImageResolution -> 300]`.

## Antes de entregar
1. Revisar balance de corchetes y `;` al final de las definiciones.
2. Revisar que cada símbolo usado esté definido antes (orden de celdas).
3. Si el código cambia el modelo, reproducir la Fig. 2A del paper (en Python
   si hace falta) y registrar en `docs/CHANGELOG.md` y `docs/VALIDACION.md`.
4. Marcar la validación como **pendiente** hasta que la persona usuaria lo
   corra.
