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
                        ── Ejemplo: el INE ──

  ETAPA 1  Inventariar todo lo que publica la fuente        (ya hecho)
           INE_TableDatasetGeneration.Rmd
                      │  recorre la API dataset a dataset
                      ▼
           INE.datasets.xlsx      el catálogo: 1.678 filas, una por dataset
                      │
  ETAPA 2  Buscar el tema dentro del catálogo, con ayuda de una IA
                      │  «de este Excel, ¿qué datasets hablan de turismo?»
                      ▼
           Dataset.Id <- 76136    el identificador del dataset elegido
                      │
  ETAPA 3  Descargar y validar ese dataset
           INE_DatasetSelection.Rmd, en dos mitades:
                      │
                      ├─ 3a  descarga  ──►  76136.rds + 76136_metadatos.xlsx
                      │
                      └─ 3b  análisis  ◄──  lee del disco esos dos ficheros
                                            y valida los requisitos AEDV
                      │
  ETAPA 4  Incorporarlo al proyecto
                      ▼
           se copian esos dos ficheros junto a la memoria y se pega en ella
           el código de la mitad 3b («Análisis del dataset»)
```

| Etapa | Qué se hace | Con qué | Qué produce | ¿La ejecuta el alumno? |
|---|---|---|---|---|
| **1. Inventario** | recorrer la API del organismo y anotar la ficha básica de **todos** sus datasets | `<FUENTE>_TableDatasetGeneration.Rmd` | el **catálogo** `<fuente>.datasets.xlsx`: una fila por dataset | **No** — ya está hecho; solo sirve para descubrir datasets nuevos |
| **2. Búsqueda del tema** | localizar en el catálogo los datasets de un tema de interés, con ayuda de una IA | el `.xlsx` del catálogo | una lista de **identificadores** candidatos | Sí |
| **3. Descarga y validación** | descargar esos datasets y comprobar si cumplen los requisitos AEDV | `<FUENTE>_DatasetSelection.Rmd` | `<nombre>.rds` (datos) + `<nombre>_metadatos.xlsx` (ficha) + el análisis de validación | Sí |
| **4. Incorporación** | llevar los datos **y el análisis** a la memoria del proyecto | los ficheros de la etapa 3 + copiar su código | la sección *Comprensión de los datos* de la memoria | Sí |

Dos detalles explican por qué el flujo está partido así:

- **El catálogo no contiene datos, solo fichas** (identificador, nombre, frecuencia,
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
