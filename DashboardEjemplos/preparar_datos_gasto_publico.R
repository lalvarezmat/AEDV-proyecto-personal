# ==============================================================================
# PREPARACIÓN DE LOS DATOS DEL CUADRO DE MANDOS Dashboard_eurostat_gasto_publico.Rmd
#
# El cuadro de mandos no lee los datos originales, sino un .rds ya preparado y
# organizado como espera la plantilla (ver el PASO 2 del cuadro de mandos): una
# tabla TIDY, con una fila por área, combinación de atributos, año e indicador,
# y las columnas
#   codigo, area, <una columna por atributo>, anio, indicador, valor
# Este script es el ejemplo de esa preparación: se ejecuta UNA vez, fuera del
# cuadro de mandos, y deja el .rds listo para leerlo con readRDS().
#
# Datos de origen: gasto de las administraciones públicas por función (COFOG)
# de Eurostat, con una columna por indicador. Metadatos en
# https://ctim.es/AEDV/data/eurostat_GASTO_PUBLICO_metadatos.xlsx
# ==============================================================================

library(tidyverse)

# código ISO de tres letras de la geometría de cada código de país de Eurostat (con EL para Grecia)
ISO3 <- c(AT = "AUT", BE = "BEL", BG = "BGR", CH = "CHE", CY = "CYP", CZ = "CZE",
          DE = "DEU", DK = "DNK", EE = "EST", EL = "GRC", ES = "ESP", FI = "FIN",
          FR = "FRA", HR = "HRV", HU = "HUN", IE = "IRL", IS = "ISL", IT = "ITA",
          LT = "LTU", LU = "LUX", LV = "LVA", MT = "MLT", NL = "NLD", NO = "NOR",
          PL = "POL", PT = "PRT", RO = "ROU", SE = "SWE", SI = "SVN", SK = "SVK")

gasto <- readRDS(url("https://ctim.es/AEDV/data/eurostat_GASTO_PUBLICO_sector.rds")) |>
  as_tibble() |>
  filter(year(time) >= 1995) |> # de 1990 a 1994 solo hay un país
  transmute(codigo = coalesce(unname(ISO3[code]), code), # código de la geometría (la UE-27 conserva el suyo)
            area = case_when(code == "ES" ~ "España", # nombre de España en español
                             code == "EU27_2020" ~ "UE-27", # nombre corto del agregado europeo
                             .default = geo), # el resto de nombres no cambia
            sector, # el atributo (con el Total como una categoría más)
            anio = as.integer(year(time)), # el año, número entero
            gasto_per_capita, pct_pib, pct_gasto, gasto) |> # indicadores (sin population, igual en todos los sectores)
  pivot_longer(c(gasto_per_capita, pct_pib, pct_gasto, gasto), # todos los indicadores en una columna...
               names_to = "indicador", values_to = "valor") |> # ...y sus valores en otra
  filter(!is.na(valor)) # pct_gasto no existe en el Total

gasto

saveRDS(gasto, "eurostat_GASTO_PUBLICO_dashboard.rds")

# GEOMETRÍA: los polígonos de los países, con las columnas codigo y nombre,
# recortados a los países de los datos (es lo que viaja al navegador en cada
# mapa). Se descarga primero a un fichero temporal porque, en algunos sistemas,
# st_read no puede acceder directamente a URLs https.
library(sf)
tmp <- tempfile(fileext = ".geojson")
curl::curl_download("https://ctim.es/AEDV/data/geo_countries.geojson", tmp)
geometria <- st_read(tmp, quiet = TRUE) |>
  select(codigo = ISO_A3, nombre = ADMIN) |>
  filter(codigo %in% gasto$codigo)
unlink(tmp)

geometria

saveRDS(geometria, "geo_gasto_publico.rds")
