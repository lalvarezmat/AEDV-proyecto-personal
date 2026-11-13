# ==============================================================================
# PREPARACIÓN DE LOS DATOS DEL CUADRO DE MANDOS Dashboard_epa_ocupados.Rmd
#
# El cuadro de mandos no lee los datos originales, sino un .rds ya preparado y
# organizado como espera la plantilla (ver el PASO 2 del cuadro de mandos): una
# tabla TIDY, con una fila por área, combinación de atributos, año e indicador,
# y las columnas
#   codigo, area, <una columna por atributo>, anio, indicador, valor
# Aquí hay DOS atributos, el sexo y el grupo de edad, cada uno en su columna y
# con su total como una categoría más («Ambos sexos» y «Total»), y DOS
# indicadores: los ocupados y los ocupados por cada 100 personas en edad de
# trabajar (16 y más años), es decir, la tasa de empleo.
# Este script es el ejemplo de esa preparación: se ejecuta UNA vez, fuera del
# cuadro de mandos, y deja los .rds listos para leerlos con readRDS().
#
# Datos de origen, de la Encuesta de Población Activa (EPA) del INE, en miles
# de personas, por comunidad autónoma, sexo y grupo de edad, con un dato por
# trimestre desde 2002 (columnas: Periodo ("2026T2"), Edad,
# Comunidades.y.Ciudades.Autónomas ("05 Canarias"), Sexo y value):
#   - EPA_65302.rds: ocupados (tabla 65302)
#   - EPA_65285.rds: población (tabla 65285). Sus grupos de edad no son
#     exactamente los de los ocupados: separa «De 55 a 64 años» y «65 y más
#     años», y su «Total» incluye a los menores de 16 años.
# ==============================================================================

library(tidyverse)

# lectura de una tabla de la EPA: los factores, a texto, y el año del trimestre
leer_epa <- function(fichero) {
  readRDS(fichero) |>
    as_tibble() |>
    mutate(across(where(is.factor), as.character)) |>
    transmute(ccaa = Comunidades.y.Ciudades.Autónomas,
              sexo = Sexo, # primer atributo
              edad = Edad, # segundo atributo
              periodo = Periodo, # el trimestre
              anio = as.integer(substr(Periodo, 1, 4)), # el año del trimestre
              value)
}

# OCUPADOS: media anual de los trimestres, solo en los años con los cuatro
ocupados <- leer_epa("EPA_65302.rds") |>
  group_by(ccaa, sexo, edad, anio) |>
  filter(n() == 4) |>
  summarise(ocupados = mean(value), .groups = "drop")

# POBLACIÓN DE 16 Y MÁS AÑOS (en edad de trabajar), con los grupos de edad de
# los ocupados: se juntan 55-64 y 65 y más, se quitan los menores de 16 y el
# Total se calcula como la suma de los grupos de 16 y más años
poblacion <- leer_epa("EPA_65285.rds") |>
  filter(!edad %in% c("Total", "Menores de 16")) |> # el Total de la tabla incluye a los menores de 16
  mutate(edad = if_else(edad %in% c("De 55 a 64 años", "65 y más años"), "De 55 y más años", edad)) |>
  group_by(ccaa, sexo, edad, periodo, anio) |>
  summarise(poblacion = sum(value), .groups = "drop")
poblacion <- bind_rows(poblacion,
                       poblacion |> # total de 16 y más años
                         group_by(ccaa, sexo, periodo, anio) |>
                         summarise(poblacion = sum(poblacion), .groups = "drop") |>
                         mutate(edad = "Total")) |>
  group_by(ccaa, sexo, edad, anio) |>
  filter(n() == 4) |> # media anual, solo en los años con los cuatro trimestres
  summarise(poblacion = mean(poblacion), .groups = "drop")

# TABLA DEL CUADRO DE MANDOS: los dos indicadores, en formato tidy
epa <- ocupados |>
  left_join(poblacion, by = c("ccaa", "sexo", "edad", "anio")) |>
  mutate(ocupados_por_100 = 100 * ocupados / poblacion) |> # ocupados por cada 100 personas de 16 y más años
  mutate(codigo = if_else(ccaa == "Total Nacional", "ES", substr(ccaa, 1, 2)), # código de la geometría ("05")
         area = if_else(ccaa == "Total Nacional", "España", sub("^[0-9]+ ", "", ccaa))) |> # nombre sin el código
  select(codigo, area, sexo, edad, anio, ocupados, ocupados_por_100) |>
  pivot_longer(c(ocupados, ocupados_por_100), names_to = "indicador", values_to = "valor") |> # tidy
  filter(!is.na(valor))

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
