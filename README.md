# Proyecto Personal — AEDV

Material del **proyecto personal** de la asignatura *Análisis Exploratorio de Datos y
Visualización* (AEDV), Grado en Ciencias de Datos, Universidad de Las Palmas de Gran Canaria.
Autor: Luis Álvarez.

**Sitio web: <https://lalvarezmat.github.io/AEDV-proyecto-personal/>**

## Qué hay aquí

| Fichero / carpeta | Para qué sirve |
|---|---|
| `README.html` | **Empieza aquí** — guía de trabajo: entregas, seguimiento y rúbrica |
| `ModeloMemoriaProyectoPersonalAEDV.Rmd` / `.html` | Plantilla de memoria que copia y rellena cada alumno |
| `RubricaEvaluaciónProyectoPersonal.xlsx` | Criterios de evaluación (27 criterios en 7 áreas) |
| `estilos.css` | Estilos de la memoria; debe acompañar al `.Rmd` |
| `SelecciónDatasetsProyecto/` | Catálogos y validadores de datasets por fuente |

En `SelecciónDatasetsProyecto/` hay una carpeta por fuente (INE, ISTAC, EUROSTAT,
WORLD_BANK, EDGAR_GHG, CLIMA_CANARIAS, CLIMA_SPAIN). Cada una contiene el documento que
genera el catálogo (`*_TableDatasetGeneration`), el que selecciona y valida el dataset
elegido (`*_DatasetSelection`) y los datos correspondientes (`.xlsx` de catálogo o `.rds`).

## Cómo obtener los ficheros

- Descarga en ZIP: <https://github.com/lalvarezmat/AEDV-proyecto-personal/archive/refs/heads/main.zip>
- O bien: `git clone https://github.com/lalvarezmat/AEDV-proyecto-personal.git`

Los `.Rmd` se abren con RStudio. Los documentos de selección de datasets leen los `.rds` y
`.xlsx` por **ruta relativa**, así que hay que conservar la estructura de carpetas.

> Aviso: en INE, ISTAC, Eurostat y Banco Mundial, la sección *Análisis del dataset* de los
> `*_DatasetSelection.Rmd` necesita los ficheros que genera la primera parte del propio
> documento; hay que hacer un *Knit* completo una vez antes de poder ejecutarla suelta.
> En EDGAR y en los dos CLIMA los datos ya vienen incluidos y compila directamente.
