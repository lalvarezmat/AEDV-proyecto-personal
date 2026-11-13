# ==============================================================================
# PREPARACIÓN DE LOS DATOS DEL CUADRO DE MANDOS Dashboard_epa_ocupados.Rmd
#
# El cuadro de mandos no lee los datos originales, sino un .rds ya preparado y
# organizado como espera la plantilla (ver el PASO 2 del cuadro de mandos): una
# tabla TIDY, con una fila por área, combinación de atributos, año e indicador,
# y las columnas
#   codigo, area, <una columna por atributo>, anio, indicador, valor
# Aquí hay DOS atributos, el sexo y el grupo de edad, cada uno en su columna y
# con su total como una categoría más («Ambos sexos» y «Total»).
# Este script es el ejemplo de esa preparación: se ejecuta UNA vez, fuera del
# cuadro de mandos, y deja los .rds listos para leerlos con readRDS().
#
# Datos de origen: EPA_65302.rds, la tabla 65302 del INE (Encuesta de Población
# Activa): ocupados, en miles de personas, por comunidad autónoma, sexo y grupo
# de edad, con un dato por trimestre desde 2002. Columnas: Periodo ("2026T2"),
# Edad, Comunidades.y.Ciudades.Autónomas ("05 Canarias"), Sexo y value.
# ==============================================================================

library(tidyverse)

epa <- readRDS("EPA_65302.rds") |>
  as_tibble() |>
  mutate(across(where(is.factor), as.character)) |> # los factores, a texto
  transmute(ccaa = Comunidades.y.Ciudades.Autónomas,
            sexo = Sexo, # primer atributo
            edad = Edad, # segundo atributo
            anio = as.integer(substr(Periodo, 1, 4)), # el año del trimestre
            value) |>
  group_by(ccaa, sexo, edad, anio) |>
  filter(n() == 4) |> # solo los años con los cuatro trimestres
  summarise(ocupados = mean(value), .groups = "drop") |> # media anual de los trimestres
  mutate(codigo = if_else(ccaa == "Total Nacional", "ES", substr(ccaa, 1, 2)), # código de la geometría ("05")
         area = if_else(ccaa == "Total Nacional", "España", sub("^[0-9]+ ", "", ccaa))) |> # nombre sin el código
  select(codigo, area, sexo, edad, anio, ocupados) |>
  pivot_longer(ocupados, names_to = "indicador", values_to = "valor") # tidy: el indicador en una columna, su valor en otra

epa

saveRDS(epa, "epa_ocupados_dashboard.rds")

# GEOMETRÍA: los polígonos de las comunidades autónomas, con las columnas
# codigo (el código INE de la comunidad, "05") y nombre. Se descarga primero a
# un fichero temporal porque, en algunos sistemas, st_read no puede acceder
# directamente a URLs https.
library(sf)
tmp <- tempfile(fileext = ".geojson")
curl::curl_download("https://ctim.es/AEDV/data/ccaa.geojson", tmp)
geometria <- st_read(tmp, quiet = TRUE) |>
  select(codigo = codauto, nombre = ccaa.shortname.es) |>
  filter(codigo %in% epa$codigo)
unlink(tmp)

geometria

saveRDS(geometria, "geo_ccaa_epa.rds")
