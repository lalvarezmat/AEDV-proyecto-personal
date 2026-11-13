# ==============================================================================
# PREPARACIÓN DE LOS DATOS DEL CUADRO DE MANDOS Dashboard_epa_ocupados.Rmd
#
# El cuadro de mandos no lee los datos originales, sino un .rds ya preparado y
# organizado como espera la plantilla (ver el PASO 2 del cuadro de mandos): una
# tabla TIDY, con una fila por área, combinación de atributos, año e indicador,
# y las columnas
#   codigo, area, <una columna por atributo>, anio, indicador, valor
# Aquí hay TRES atributos, el sexo, el grupo de edad y el sector económico,
# cada uno en su columna y con su total como una categoría más («Ambos sexos»,
# «Total» y «Total»), y DOS indicadores: los ocupados y los ocupados por cada
# 100 personas en edad de trabajar (16 y más años) del mismo sexo y grupo de
# edad (con el sector Total, la tasa de empleo).
# Este script es el ejemplo de esa preparación: se ejecuta UNA vez, fuera del
# cuadro de mandos, y deja los .rds listos para leerlos con readRDS().
#
# Datos de origen, de la Encuesta de Población Activa (EPA) del INE, en miles
# de personas, por comunidad autónoma, con un dato por trimestre:
#   - EPA_65304.rds: ocupados por grupo de edad, sexo y sector económico
#     (tabla 65304, desde 2008). El INE no publica el dato cuando la muestra es
#     pequeña (por ejemplo, los jóvenes en la agricultura): esas celdas son NA.
#   - EPA_65285.rds: población por grupo de edad y sexo (tabla 65285). Sus
#     grupos de edad son más finos que los de los ocupados (25-34, 35-44,
#     45-54; 55-64 y 65 y más) y su «Total» incluye a los menores de 16 años.
# ==============================================================================

library(tidyverse)
library(sf)

# GEOMETRÍA: los polígonos de las comunidades autónomas, con las columnas
# codigo (el código INE de la comunidad, "05") y nombre (el nombre corto:
# "Madrid", "Baleares"). Se lee la primera porque sus nombres cortos son
# también los de las áreas de los datos. Se descarga primero a un fichero
# temporal porque, en algunos sistemas, st_read no puede acceder directamente
# a URLs https.
tmp <- tempfile(fileext = ".geojson")
curl::curl_download("https://ctim.es/AEDV/data/ccaa.geojson", tmp)
geometria <- st_read(tmp, quiet = TRUE) |>
  select(codigo = codauto, nombre = ccaa.shortname.es)
unlink(tmp)

# nombre corto de cada comunidad a partir de su código
NOMBRES <- setNames(geometria$nombre, geometria$codigo)

# lectura de una tabla de la EPA: los factores, a texto, y el año del trimestre
leer_epa <- function(fichero) {
  readRDS(fichero) |>
    as_tibble() |>
    mutate(across(where(is.factor), as.character)) |>
    rename(ccaa = Comunidades.y.Ciudades.Autónomas, sexo = Sexo, edad = Edad) |>
    mutate(anio = as.integer(substr(Periodo, 1, 4))) # el año del trimestre
}

# OCUPADOS: media anual de los trimestres, solo con los cuatro trimestres
# publicados (sin NA)
ocupados <- leer_epa("EPA_65304.rds") |>
  rename(sector = Sector.económico) |>
  filter(!is.na(value)) |>
  group_by(ccaa, sexo, edad, sector, anio) |>
  filter(n() == 4) |>
  summarise(ocupados = mean(value), .groups = "drop")

# POBLACIÓN DE 16 Y MÁS AÑOS (en edad de trabajar), con los grupos de edad de
# los ocupados: 25-34, 35-44 y 45-54 se juntan en 25-54; 55-64 y 65 y más, en
# 55 y más; se quitan los menores de 16, y el Total se calcula como la suma de
# los grupos de 16 y más años
poblacion <- leer_epa("EPA_65285.rds") |>
  filter(!edad %in% c("Total", "Menores de 16")) |> # el Total de la tabla incluye a los menores de 16
  mutate(edad = case_when(edad %in% c("De 25 a 34 años", "De 35 a 44 años", "De 45 a 54 años") ~ "De 25 a 54 años",
                          edad %in% c("De 55 a 64 años", "65 y más años") ~ "De 55 y más años",
                          .default = edad)) |>
  group_by(ccaa, sexo, edad, Periodo, anio) |>
  summarise(poblacion = sum(value), .groups = "drop")
poblacion <- bind_rows(poblacion,
                       poblacion |> # total de 16 y más años
                         group_by(ccaa, sexo, Periodo, anio) |>
                         summarise(poblacion = sum(poblacion), .groups = "drop") |>
                         mutate(edad = "Total")) |>
  group_by(ccaa, sexo, edad, anio) |>
  filter(n() == 4) |> # media anual, solo en los años con los cuatro trimestres
  summarise(poblacion = mean(poblacion), .groups = "drop")

# TABLA DEL CUADRO DE MANDOS: los dos indicadores, en formato tidy
epa <- ocupados |>
  left_join(poblacion, by = c("ccaa", "sexo", "edad", "anio")) |>
  mutate(ocupados_por_100 = 100 * ocupados / poblacion) |> # por cada 100 personas de 16 y más años del mismo sexo y edad
  mutate(codigo = if_else(ccaa == "Total Nacional", "ES", substr(ccaa, 1, 2)), # código de la geometría ("05")
         area = if_else(codigo == "ES", "España", unname(NOMBRES[codigo]))) |> # nombre corto ("Madrid", no "13 Madrid, Comunidad de")
  select(codigo, area, sexo, edad, sector, anio, ocupados, ocupados_por_100) |>
  pivot_longer(c(ocupados, ocupados_por_100), names_to = "indicador", values_to = "valor") |> # tidy
  filter(!is.na(valor)) |>
  # el cuadro de mandos ofrece las categorías en el orden en que aparecen: el del INE
  arrange(match(sector, c("Agricultura", "Industria", "Construcción", "Servicios", "Total")),
          match(edad, c("De 16 a 19 años", "De 20 a 24 años", "De 25 a 54 años", "De 55 y más años", "Total")),
          match(sexo, c("Hombres", "Mujeres", "Ambos sexos")),
          codigo, anio, indicador)

epa

saveRDS(epa, "epa_ocupados_dashboard.rds")

# la geometría, recortada a las comunidades de los datos
geometria <- geometria |>
  filter(codigo %in% epa$codigo)

geometria

saveRDS(geometria, "geo_ccaa_epa.rds")
