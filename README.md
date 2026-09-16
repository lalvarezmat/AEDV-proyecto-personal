# Proyecto Personal — AEDV

Material del **proyecto personal** de la asignatura *Análisis Exploratorio de Datos y
Visualización* (AEDV), Grado en Ciencia e Ingeniería de Datos, Universidad de Las Palmas de Gran Canaria.
Autor: Luis Álvarez.

**Sitio web: <https://lalvarezmat.github.io/AEDV-proyecto-personal/>**

## Qué hay aquí

| Fichero / carpeta | Para qué sirve |
|---|---|
| `README_Proyecto_Personal.html` | **Empieza aquí** — guía de trabajo: entregas, seguimiento y rúbrica |
| `ModeloMemoriaProyectoPersonalAEDV.Rmd` / `.html` | Plantilla de memoria que copia y rellena cada alumno |
| `RubricaEvaluaciónProyectoPersonal.xlsx` | Criterios de evaluación (27 criterios en 7 áreas) |
| `estilos.css` | Estilos de la memoria; debe acompañar al `.Rmd` |
| `SelecciónDatasetsProyecto/` | Catálogos y validadores de datasets por fuente |

En `SelecciónDatasetsProyecto/` hay una carpeta por fuente (INE, ISTAC, EUROSTAT,
WORLD_BANK, EDGAR_GHG, CLIMA_CANARIAS, CLIMA_SPAIN, INFORME_PISA). Cada una contiene el
documento que genera el catálogo (`*_TableDatasetGeneration`), el que descarga, valida y
analiza los datasets elegidos (`*_DatasetSelection`) y los datos correspondientes
(`.xlsx` de catálogo o `.rds`).

## Cómo se eligen los datasets: las cuatro etapas

Las cuatro fuentes con API —**INE**, **ISTAC**, **EUROSTAT** y **Banco Mundial**— se
abordan **exactamente igual**. Cambia el organismo y el nombre de los ficheros, no el
procedimiento; lo que produce cada etapa es justo de lo que parte la siguiente:

```
  ETAPA 1   catálogo (inventario) de todos los datasets, que se guarda en la hoja excel <FUENTE>.datasets.xlsx (ya está hecho)
                     ▼
  ETAPA 2   buscas en ese Excel, con ayuda de una IA, los datasets de tu tema
                     ▼
  ETAPA 3   descarga y valida los que has elegido usando <FUENTE>_DatasetSelection.Rmd
                     ▼
  ETAPA 4   copias a tu memoria de proyecto personal los ficheros que ha guardado y su análisis
```

| Etapa | Con qué se hace y qué produce | ¿La ejecuta el alumno? |
|---|---|---|
| **1. Inventario** *(ya hecho)* — recorrer la API del organismo y anotar la ficha básica de **todos** sus datasets | `<FUENTE>_TableDatasetGeneration.Rmd` → el **catálogo** `<FUENTE>.datasets.xlsx`, una fila por dataset | No |
| **2. Búsqueda del tema** — localizar en el catálogo los datasets de un tema de interés | el `.xlsx` del catálogo y una IA → una lista de **identificadores** candidatos | Sí |
| **3. Descarga y validación** — descargar esos datasets y comprobar si cumplen los requisitos AEDV | `<FUENTE>_DatasetSelection.Rmd` → `<nombre>.rds` (datos) y `<nombre>_metadatos.xlsx` (ficha), más el análisis de validación | Sí |
| **4. Incorporación** — llevar los datos **y el análisis** a la memoria del proyecto | los ficheros de la etapa 3 y su código → la sección *Comprensión de los datos* de la memoria | Sí |

Dos detalles explican por qué el flujo está partido así:

- **El catálogo, con el inventario, no contiene datos, solo fichas** (identificador, nombre, frecuencia,
  cobertura geográfica, fechas, dimensiones). Los datos no se descargan hasta la etapa 3.
- **La etapa 3 descarga antes de analizar, a propósito.** El `*_DatasetSelection.Rmd`
  primero guarda en disco el `.rds` y el `_metadatos.xlsx` (3a) y después analiza
  **leyendo esos ficheros** (3b). Por eso el código de la sección *Análisis del dataset*
  se puede pegar tal cual en la memoria —arranca leyendo del disco— y por eso el proyecto
  deja de depender de que el organismo mantenga el dataset disponible.

EDGAR, CLIMA e INFORME PISA son datasets únicos ya construidos: **no tienen etapas 1 y 2**
(no hay catálogo en el que buscar), pero conservan la misma estructura 3a/3b.

El detalle completo está en
[`SelecciónDatasetsProyecto/README_Selección_Tema.html`](https://lalvarezmat.github.io/AEDV-proyecto-personal/Selecci%C3%B3nDatasetsProyecto/README_Selecci%C3%B3n_Tema.html).

## Cómo obtener los ficheros

- Descarga en ZIP: <https://github.com/lalvarezmat/AEDV-proyecto-personal/archive/refs/heads/main.zip>
- O bien: `git clone https://github.com/lalvarezmat/AEDV-proyecto-personal.git`

Los `.Rmd` se abren con RStudio. Los documentos de selección de datasets leen los `.rds` y
`.xlsx` por **ruta relativa**, así que hay que conservar la estructura de carpetas.

> Aviso: en INE, ISTAC, Eurostat y Banco Mundial, la sección *Análisis del dataset* de los
> `*_DatasetSelection.Rmd` necesita los ficheros que genera la primera parte del propio
> documento; hay que hacer un *Knit* completo una vez antes de poder ejecutarla suelta.
> En EDGAR y en los dos CLIMA los datos ya vienen incluidos y compila directamente.
