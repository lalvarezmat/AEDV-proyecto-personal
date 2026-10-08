# ==========================================================================
# Funciones auxiliares del curso AEDV
# ==========================================================================
#
# Este es el único fichero de utilidades del curso: se alimentan de él tanto
# los capítulos del libro como los cuadros de mando de ejemplo. Los documentos
# lo cargan con source(), bien desde una ruta relativa dentro del proyecto
# (source("../utilidades.R")), bien desde el servidor del curso
# (source("https://ctim.es/AEDV/utilidades.R")) cuando el fichero se distribuye
# suelto al alumno, como ocurre con los dashboards.
#
# Todas las funciones propias del curso llevan el prefijo aedv_ para que se
# distingan a simple vista de las de librería. Las que empiezan por un punto
# (.aedv_) son auxiliares internas: las usan otras funciones de este fichero y
# no están pensadas para llamarse directamente.
#
# Las funciones asumen que tidyverse está cargado en la sesión (usan tibble,
# dplyr y el pipe nativo), lo cual se cumple en todos los documentos del curso.

library(utils)    # adist(): distancia de edición entre cadenas de texto
library(leaflet)  # mapas interactivos, para el mapa coroplético del final


# ==========================================================================
# Cálculo numérico y transformaciones
# ==========================================================================

# Función que, dado un vector de números y un nivel de significancia de Pareto, calcula el valor
# del vector tal que la suma de los valores por encima de ese valor supera el total de la suma
# multiplicado por el nivel de significancia.
aedv_pareto_value <- function(
v, # vector numérico
ParetoSignificancia # número en [0,1] con el grado de significancia de Pareto
){

  x <- sort(v,decreasing = TRUE)
  return(x[which(cumsum(x)>=ParetoSignificancia*sum(x))[1]])
}


# Calcula la derivada de una serie temporal
aedv_derivada <-function(data){
  res <- numeric(length(data))
  res[1] <- NA
  for(i in 2:length(res) ){
    res[i]=data[i]-data[i-1]
  }
  return(res)
}


# Cálculo de la transformación de Yeo–Johnson a partir de un vector y un valor de lambda.
# Los valores ausentes se conservan como tales (antes se convertían en ceros).
aedv_yeo_johnson <- function(y, lambda) {
  y_t <- rep(NA_real_, length(y))

  # Para y >= 0
  pos_idx <- which(y >= 0)
  if (lambda == 0) {
    y_t[pos_idx] <- log(y[pos_idx] + 1)
  } else {
    y_t[pos_idx] <- ((y[pos_idx] + 1)^lambda - 1) / lambda
  }

  # Para y < 0
  neg_idx <- which(y < 0)
  if (lambda == 2) {
    y_t[neg_idx] <- -log(-y[neg_idx] + 1)
  } else {
    y_t[neg_idx] <- -(((-y[neg_idx] + 1)^(2 - lambda) - 1) / (2 - lambda))
  }

  return(y_t)
}


# Estimación óptima de lambda para la transformación de Yeo–Johnson a un vector y para el R2 de la regresión lineal con un vector x sea máximo
aedv_lambda_optimo <- function(x, y, lambda_range = c(-1, 1.9)) {

  # Función objetivo: R² negativo (porque optimize minimiza)
  r2_neg <- function(lambda) {
    y_t <- aedv_yeo_johnson(y, lambda)
    modelo <- lm(y_t ~ x)
    return(-summary(modelo)$r.squared)  # queremos maximizar R²
  }

  # Optimización de lambda
  opt <- optimize(r2_neg, interval = lambda_range)

  # Se retorna el valor óptimo de lambda
  return(opt$minimum)
}


# Estimación de lambda para la transformación de Yeo-Johnson de un vector x que
# MAXIMIZA LA VEROSIMILITUD de que los valores transformados sigan una normal,
# es decir, el lambda que mejor aproxima x a una distribución normal.
# A diferencia de aedv_lambda_optimo(), que busca el mejor ajuste lineal con
# OTRA variable, aquí solo interviene x. La búsqueda se acota a un intervalo
# razonable: para una variable ya casi simétrica el óptimo sin acotar puede
# dispararse a valores extremos sin mejorar apenas la asimetría.
aedv_yeojohnson_lambda <- function(x, intervalo = c(-1, 3)) {
  x <- as.numeric(x[!is.na(x)])

  # Log-verosimilitud del modelo de Yeo-Johnson para un valor de lambda
  log_verosimilitud <- function(lambda) {
    z <- aedv_yeo_johnson(x, lambda)
    n <- length(z)
    -n/2 * log(mean((z - mean(z))^2)) + (lambda - 1) * sum(sign(x) * log1p(abs(x)))
  }

  optimize(log_verosimilitud, interval = intervalo, maximum = TRUE)$maximum
}


# Ajusta la recta de regresión de y sobre x y devuelve, en una tabla de una
# fila, los indicadores habituales para valorar la calidad del ajuste:
#   sigma         : desviación típica de los residuos (coincide con el RMSE)
#   R2            : coeficiente de determinación
#   F_p_value     : p-valor del estadístico F (significancia global del modelo)
#   SW_p_value    : p-valor del test de normalidad de Shapiro-Wilk sobre los
#                   residuos (valores pequeños indican que no son normales)
#   BP_p_value    : p-valor del test de Breusch-Pagan (valores pequeños indican
#                   heterocedasticidad: la dispersión de los residuos no es
#                   constante). Se calcula como n·R² de la regresión auxiliar de
#                   los residuos al cuadrado sobre x, que es exactamente lo que
#                   devuelve lmtest::bptest(), sin necesidad de esa librería
#   pendiente, termino_indep : coeficientes de la recta y = a·x + b
aedv_lm_analisis <- function(x, y) {
  datos <- data.frame(x = as.numeric(x), y = as.numeric(y))
  datos <- datos[stats::complete.cases(datos), ]

  modelo <- lm(y ~ x, data = datos)
  resumen <- summary(modelo)
  f <- resumen$fstatistic

  # Test de Breusch-Pagan: regresión auxiliar de los residuos al cuadrado
  aux <- lm(residuals(modelo)^2 ~ datos$x)
  bp <- nrow(datos) * summary(aux)$r.squared

  data.frame(
    sigma         = resumen$sigma,
    R2            = resumen$r.squared,
    F_p_value     = pf(f[1], f[2], f[3], lower.tail = FALSE),
    SW_p_value    = shapiro.test(residuals(modelo))$p.value,
    BP_p_value    = pchisq(bp, df = 1, lower.tail = FALSE),
    pendiente     = coef(modelo)[2],
    termino_indep = coef(modelo)[1],
    row.names     = NULL
  )
}


# ==========================================================================
# Diagnóstico de residuos
# ==========================================================================
#
# Un modelo bien ajustado deja unos residuos que se comportan como RUIDO
# BLANCO: media cero, sin autocorrelación, normales y de varianza constante.
# Comprobarlo supone lanzar cuatro contrastes distintos y leer cuatro salidas
# con formatos diferentes, lo que hace difícil ver el diagnóstico de conjunto.
# Estas dos funciones lo reducen a dos llamadas: aedv_tests_residuos() calcula
# los cuatro y devuelve una tabla, y aedv_tabla_tests() la presenta coloreando
# el p-valor en verde o rojo según se supere o no cada contraste.

# Calcula de una vez los cuatro tests de diagnóstico sobre un vector de
# residuos. Devuelve un tibble con una fila por test y las columnas:
#   propiedad : la propiedad del ruido blanco que se está contrastando
#   test      : nombre del contraste utilizado
#   p_valor   : su p-valor (NA si el test no es aplicable, ver más abajo)
#   supera    : TRUE si p_valor > alpha, es decir, si NO hay evidencia para
#               rechazar la propiedad
#
# Argumentos:
#   residuos : vector numérico con los residuos (o el remainder de una
#              descomposición). Los NA se descartan.
#   tiempo   : vector con el instante de cada residuo, necesario para el test
#              de Breusch-Pagan; si se omite se usa el orden de la serie.
#   lag      : nº de retardos del test de Ljung-Box.
#   alpha    : nivel de significación con el que se decide "supera".
#   dof      : grados de libertad que resta Ljung-Box (nº de parámetros del
#              modelo). Se deja en 0 por omisión, igual que feasts::ljung_box.
#
# NOTA sobre Shapiro-Wilk: stats::shapiro.test() solo admite entre 3 y 5000
# observaciones. Fuera de ese rango se devuelve NA en lugar de interrumpir el
# diagnóstico, porque con muestras grandes ese test deja de ser informativo
# (rechaza la normalidad ante desviaciones minúsculas) y conviene sustituirlo
# por el diagnóstico gráfico del QQ-plot.
aedv_tests_residuos <- function(residuos, tiempo = NULL, lag = 10,
                                alpha = 0.05, dof = 0) {
  residuos <- as.numeric(residuos)
  if (is.null(tiempo)) tiempo <- seq_along(residuos)

  # los NA se descartan de forma conjunta para que residuos y tiempo sigan
  # alineados (los residuos de un modelo ajustado suelen tener NA al principio)
  ok <- !is.na(residuos)
  residuos <- residuos[ok]
  tiempo   <- as.numeric(tiempo)[ok]
  n <- length(residuos)

  p_sw <- if (n >= 3 && n <= 5000) stats::shapiro.test(residuos)$p.value else NA_real_

  tibble(
    propiedad = c("Media cero", "Ausencia de autocorrelación",
                  "Normalidad", "Homocedasticidad (varianza constante)"),
    test      = c("t de Student", "Ljung-Box", "Shapiro-Wilk", "Breusch-Pagan"),
    p_valor   = c(
      stats::t.test(residuos)$p.value,
      stats::Box.test(residuos, lag = lag, type = "Ljung-Box", fitdf = dof)$p.value,
      p_sw,
      lmtest::bptest(residuos ~ tiempo)$p.value
    )
  ) |>
    mutate(supera = .data$p_valor > alpha)
}


# Formatea un p-valor para la tabla: notación científica cuando es diminuto,
# cuatro decimales en caso contrario, y un texto explícito cuando no aplica.
.aedv_formato_p <- function(p) {
  ifelse(is.na(p), "no aplicable",
         ifelse(p < 1e-4, formatC(p, format = "e", digits = 1),
                formatC(p, format = "f", digits = 4)))
}


# Presenta el resultado de aedv_tests_residuos() como una tabla gt en la que el
# p-valor aparece sobre fondo VERDE si el test se supera y ROJO si no, de modo
# que el diagnóstico completo se lee de un vistazo. Los tests no aplicables
# (Shapiro-Wilk con más de 5000 datos) quedan en gris.
aedv_tabla_tests <- function(tests, titulo = "Diagnóstico de residuos",
                             alpha = 0.05) {
  d <- tests |>
    mutate(
      p_txt     = .aedv_formato_p(.data$p_valor),
      cond      = paste0("p-valor > ", alpha),
      veredicto = case_when(is.na(.data$supera) ~ "—",
                            .data$supera        ~ "Sí",
                            TRUE                ~ "No")) |>
    select("propiedad", "test", "cond", "p_txt", "veredicto")

  filas_ok <- which(tests$supera %in% TRUE)   # verde
  filas_ko <- which(tests$supera %in% FALSE)  # rojo
  filas_na <- which(is.na(tests$supera))      # gris: test no aplicable

  # id HTML fijo derivado del título: el aleatorio de gt se repetía entre tablas
  # (mismo estado del generador aleatorio) y dejaba ids duplicados en la página
  id_tabla <- paste0("tests-", gsub("[^a-z0-9]+", "-",
                     tolower(iconv(titulo, to = "ASCII//TRANSLIT"))))

  tb <- d |>
    gt::gt(id = id_tabla) |>
    gt::tab_header(title = titulo) |>
    gt::cols_label(propiedad = "Propiedad que se analiza",
                   test      = "Test",
                   cond      = "Condición para superarlo",
                   p_txt     = "p-valor",
                   veredicto = "¿Lo supera?") |>
    gt::cols_align("center", columns = c("p_txt", "veredicto")) |>
    gt::tab_options(table.font.size = gt::px(13))

  pinta <- function(tb, filas, color) {
    if (length(filas) == 0) return(tb)
    tb |> gt::tab_style(
      style = list(gt::cell_fill(color = color), gt::cell_text(weight = "bold")),
      locations = gt::cells_body(columns = c("p_txt", "veredicto"), rows = filas))
  }

  tb |>
    pinta(filas_ok, "#c8e6c9") |>
    pinta(filas_ko, "#ffcdd2") |>
    pinta(filas_na, "#e0e0e0")
}


# ==========================================================================
# Emparejamiento aproximado de textos
# ==========================================================================
#
# Combinar dos tablas por un nombre de texto (municipios, países) cuando no hay
# un código identificativo común. Se hace SIEMPRE en tres pasos, porque el paso
# intermedio no se puede automatizar:
#
#   1. aedv_left_join_nearest_string() empareja los valores DISTINTOS de la
#      columna de texto y devuelve, para cada uno, el más parecido de la otra
#      tabla junto con la distancia.
#   2. Se revisa ese resultado a mano, consultando alternativas con
#      aedv_candidatos(), y se corrige lo que esté mal. Al consultar
#      candidatos para un texto dudoso, hay que pasarle el VECTOR RESTANTE
#      (los valores de la otra tabla que ningún texto ha emparejado ya de
#      forma perfecta, típicamente con setdiff()): un valor ya reservado por
#      un matching perfecto no puede ser también la respuesta correcta de
#      otro texto distinto.
#   3. aedv_left_join_completo() aplica los emparejamientos ya verificados a las
#      tablas completas.

# Parecido entre dos textos como porcentaje (0 a 100), a partir de la distancia
# de inserción/borrado normalizada por la longitud de los dos textos.
.aedv_ratio <- function(a, b) {
  if (nchar(a) + nchar(b) == 0) return(100)
  d <- adist(a, b, costs = c(insertions = 1, deletions = 1, substitutions = 2))[1, 1]
  100 * (1 - d / (nchar(a) + nchar(b)))
}

# Palabras distintas de un texto, en minúsculas y ordenadas alfabéticamente
.aedv_tokens <- function(s) {
  t <- unlist(strsplit(tolower(s), "[^[:alnum:]áéíóúüñ]+"))
  sort(unique(t[t != ""]))
}

# Parecido entre los CONJUNTOS DE PALABRAS de dos textos (0 a 100), ignorando el
# orden y las palabras sobrantes: "Rioja, La" y "La Rioja" dan 100, y también lo
# da un texto contenido en otro, como "San Miguel" dentro de "San Miguel de Abona".
.aedv_token_set_ratio <- function(a, b) {
  ta <- .aedv_tokens(a)
  tb <- .aedv_tokens(b)
  comunes <- paste(intersect(ta, tb), collapse = " ")
  # al texto común se le añade lo que sobra de cada lado, y se comparan los tres
  solo_a <- trimws(paste(comunes, paste(setdiff(ta, tb), collapse = " ")))
  solo_b <- trimws(paste(comunes, paste(setdiff(tb, ta), collapse = " ")))
  max(.aedv_ratio(comunes, solo_a),
      .aedv_ratio(comunes, solo_b),
      .aedv_ratio(solo_a, solo_b))
}


# PASO 1 del emparejamiento aproximado.
# Función que hace un "left join" entre las columnas de texto de dos tablas u
# y v usando adist() de la librería utils. Devuelve un tibble con:
#   pos1   : la posición en u
#   pos2   : la posición en v del string más cercano a u[pos1]
#   value1 : el valor de u
#   value2 : el valor de v más cercano a value1
#   dis    : la distancia entre value1 y value2 (0 = matching perfecto)
# Nota : hay que tener en cuenta que la comparación puede fallar
# y por tanto es conveniente verificar los resultados.
#
# Se hace en DOS PASADAS para reducir los errores del emparejamiento
# automático:
#   1. Se busca, para cada texto de u, el más cercano de TODO v, sin
#      restricciones. Los que salen con dis==0 (matching perfecto) reservan
#      su texto de v en exclusiva: ningún otro texto de u puede tener ese
#      mismo texto como el correcto, así que no tiene sentido ofrecérselo.
#   2. Para los que NO han tenido un matching perfecto (dis>0, los que de
#      todas formas hay que revisar a mano) se recalcula su mejor candidato,
#      pero esta vez SOLO entre los textos de v que ningún otro texto de u
#      se ha quedado ya de forma perfecta. Sin este segundo cálculo, un texto
#      sin correspondencia real (una fila agregada, un total) puede aparecer
#      "emparejado" con un texto de v que en realidad es la respuesta
#      correcta de otro texto distinto, lo cual confunde la revisión manual.
aedv_left_join_nearest_string <- function(u, v, ignore_case = FALSE) {
  textos_u <- as.character(pull(u, 1))
  textos_v <- as.character(pull(v, 1))

  m <- adist(textos_u, textos_v, ignore.case = ignore_case) # fila=u, columna=v

  # PRIMERA PASADA: el texto de v más cercano a cada texto de u, sin restringir
  mejor <- apply(m, 1, which.min)
  pos2 <- mejor
  dis  <- m[cbind(seq_along(textos_u), mejor)]

  # los matching PERFECTOS (dis==0) reservan su texto de v en exclusiva
  perfectos <- dis == 0
  usados <- unique(pos2[perfectos])

  # SEGUNDA PASADA: para los que no han tenido matching perfecto, el mejor
  # candidato SOLO entre los textos de v libres (sin matching perfecto de
  # otro texto de u)
  if (length(usados) > 0 && any(!perfectos) && length(usados) < length(textos_v)) {
    libres <- setdiff(seq_along(textos_v), usados)
    m_libres <- m[!perfectos, libres, drop = FALSE]
    mejor_libre <- apply(m_libres, 1, which.min)
    pos2[!perfectos] <- libres[mejor_libre]
    dis[!perfectos]  <- m_libres[cbind(seq_len(nrow(m_libres)), mejor_libre)]
  }

  tibble(
    pos1   = seq_along(textos_u),
    pos2   = pos2,
    value1 = textos_u,
    value2 = textos_v[pos2],
    dis    = as.integer(dis)
  )
}


# PASO 2 del emparejamiento aproximado.
# Lista los n valores de v más parecidos a un texto, para decidir a mano un
# emparejamiento que aedv_left_join_nearest_string ha fallado. Devuelve un tibble
# con el candidato y DOS medidas de parecido, porque cada una falla en casos
# distintos:
#   dis       : distancia de edición (0 = idénticos), la que usa el emparejamiento
#               automático. Penaliza mucho las diferencias de longitud, así que un
#               nombre corto y ajeno puede ganar al correcto pero más largo.
#   similitud : parecido entre los conjuntos de palabras (0 a 100), que ignora el
#               orden y las palabras sobrantes.
# Conviene mirar las dos columnas: la decisión final es SIEMPRE de la persona.
# IMPORTANTE: al revisar un resultado de aedv_left_join_nearest_string(), pasar
# como v el VECTOR RESTANTE (setdiff(v_completo, res$value2[res$dis == 0])), no
# el vector completo: un valor con matching perfecto ya es la respuesta correcta
# de otro texto, así que no debe ofrecerse como candidato de este.
aedv_candidatos <- function(texto, v, n = 5, ignore_case = FALSE) {
  valores <- as.character(v[!is.na(v)])
  norm <- if (ignore_case) tolower else identity
  consulta <- norm(as.character(texto))
  normalizados <- norm(valores)

  tibble(
    pos = seq_along(valores),
    candidato = valores,
    dis = as.integer(adist(consulta, normalizados)[1, ]),
    similitud = round(vapply(normalizados,
                             function(s) .aedv_token_set_ratio(consulta, s),
                             numeric(1)), 1)
  ) |>
    arrange(desc(similitud), dis) |>
    head(n)
}


# PASO 3 del emparejamiento aproximado.
# Combina dos tablas usando una tabla de emparejamientos "res" ya VERIFICADA (y
# corregida a mano si hacía falta): cambiando value2 por el valor correcto, o
# poniéndolo a NA para declarar que ese texto NO tiene correspondencia.
# Las filas de izq cuyo texto se anuló (value2 a NA) se conservan, con las
# columnas de der a NA: es una combinación por la izquierda, no un filtro. Antes
# de combinar se comprueba que los emparejamientos no vayan a MULTIPLICAR filas,
# que es el fallo silencioso más habitual al combinar tablas.
aedv_left_join_completo <- function(izq, der, res, izq_on, der_on) {
  if (!all(c("value1", "value2") %in% names(res))) {
    stop("res debe tener las columnas 'value1' y 'value2' ",
         "(la salida de aedv_left_join_nearest_string)")
  }

  mapa <- res |>
    filter(!is.na(value2)) |>
    select(value1, value2)

  # cada texto debe aparecer una sola vez, o la combinación multiplicaría filas
  repetidos <- unique(mapa$value1[duplicated(mapa$value1)])
  if (length(repetidos) > 0) {
    stop("hay ", length(repetidos), " texto(s) repetidos en res$value1 (",
         paste(head(repetidos, 3), collapse = ", "),
         "...): cada texto debe aparecer una sola vez, o la combinación ",
         "multiplicaría filas. Empareja los valores únicos de la columna, ",
         "no todas sus filas.")
  }

  # los nombres homónimos no identifican una única fila de la tabla de destino
  valores_der <- as.character(der[[der_on]])
  ambiguos <- unique(valores_der[duplicated(valores_der) |
                                 duplicated(valores_der, fromLast = TRUE)])
  choque <- sort(intersect(as.character(mapa$value2), ambiguos))
  if (length(choque) > 0) {
    stop(length(choque), " valor(es) de der$", der_on, " están repetidos y no ",
         "identifican una única fila (", paste(head(choque, 3), collapse = ", "),
         "...). Con nombres homónimos hay que combinar por un código ",
         "identificativo, no por el texto.")
  }

  izq |>
    mutate(.texto = as.character(.data[[izq_on]])) |>
    left_join(mapa, by = c(".texto" = "value1")) |>
    left_join(der, by = c("value2" = der_on)) |>
    select(-.texto, -value2)
}


# ==========================================================================
# Visualización
# ==========================================================================

# Dibuja un gráfico de líneas con DOBLE EJE VERTICAL: la segunda variable se
# reescala al rango de la primera para poder superponerlas, y el eje de la
# derecha deshace esa escala para que sus valores se lean en su unidad real.
# Ojo: un doble eje puede sugerir relaciones que no existen, porque la escala de
# cada eje se elige libremente. Ver la discusión del capítulo de visualización
# estática antes de usarlo.
aedv_doble_eje <- function(
tb,     # tibble con los datos
namex,  # nombre de la variable del eje x
namey1, # nombre de la variable del eje y por la izquierda
color1, # color con el que se dibuja el eje y de la izquierda
namey2, # nombre de la variable del eje y por la derecha
color2  # color con el que se dibuja el eje y de la derecha
){
  y1.min <- min(tb[[namey1]])
  y2.min <- min(tb[[namey2]])
  m1 <- (max(tb[[namey1]]) - y1.min) / (max(tb[[namey2]]) - y2.min)
  p <- tibble(
    x  = tb[[namex]],
    y1 = tb[[namey1]],
    y2 = tb[[namey2]]
  ) |> ggplot(aes(x = x)) +
  geom_line(aes(y = y1, color = color1)) +                      # primera variable
  geom_line(aes(y = y1.min + (y2 - y2.min) * m1), color = color2) +  # segunda, reescalada
  # títulos de los ejes; el eje derecho deshace el reescalado
  scale_y_continuous(name = namey1,
                     sec.axis = sec_axis(~ (. - y1.min)/m1 + y2.min, name = namey2)) +
  theme(
    axis.title.y = element_text(color = color1),
    axis.title.y.right = element_text(color = color2),
    legend.position = "none"
  )
  return(p)
}


# Dibuja una MATRIZ DE CORRELACIÓN interactiva con plotly, con una escala de
# color rojo (-1) → blanco (0) → azul (+1) anclada a esos extremos. Es la
# alternativa a los gráficos estáticos de correlación cuando el resultado va a
# ir dentro de un cuadro de mando.
library(plotly)
aedv_matriz_correlacion <- function(
M,                                  # matriz de correlación
show_values = TRUE,                 # imprimir los valores dentro de las celdas
title = "matriz de correlación"     # título del gráfico
) {

  # Texto dentro de la celda (2 decimales)
  cell_text <- round(M, 2)

  # Escala personalizada: rojo -> blanco -> azul
  custom_colors <- list(
    c(0, "red"),      # mínimo (-1) → rojo
    c(0.5, "white"),  # 0 → blanco
    c(1, "steelblue") # máximo (+1) → azul
  )

  # Crear heatmap interactivo
  plot_ly(
    x = colnames(M),
    y = rownames(M),
    z = M,
    type = "heatmap",
    colorscale = custom_colors,
    zmin = -1,
    zmax = 1,
    text = if (show_values) cell_text else NULL,
    texttemplate = if (show_values) "%{text}" else NULL,
    textfont = list(color = "black"),
    hovertemplate = "X: %{x}<br>Y: %{y}<br>Correlación: %{z:.4f}<extra></extra>"
  ) |>
    layout(
      title = title,
      xaxis = list(title = "", tickangle = 45),
      yaxis = list(title = "", autorange = "reversed")
    )
}


# Lee un GeoJSON remoto y lo devuelve como objeto sf en coordenadas WGS84.
# Se descarga primero a un fichero temporal porque, en algunos sistemas Windows,
# st_read() no puede acceder directamente a URLs https.
aedv_leer_geojson <- function(url) {
  tmp <- tempfile(fileext = ".geojson")
  curl::curl_download(url, tmp)
  st_read(tmp, quiet = TRUE) |>
    st_transform(4326) # normalizamos a WGS84 (longitud/latitud)
}


# Dibuja la PREDICCIÓN de un modelo fable junto con la serie observada y sus
# intervalos de confianza al 80% y 95%. Se usa en lugar de autoplot() para
# evitar conflictos con patchwork y controlar el estilo del gráfico.
# Si se pasa el modelo ajustado en "fit", superpone además su curva de valores
# ajustados sobre la parte observada, para ver cómo reproduce la serie.
aedv_dibujar_prediccion <- function(fc, historico, fit = NULL) {
  fc_ci <- fc |>
    hilo(level = c(80, 95)) |>
    unpack_hilo(c(`80%`, `95%`))
  # Extraer nombre de la variable y del índice temporal
  y_var <- measured_vars(fc)[1]
  t_var <- as.character(tsibble::index(fc))
  # Dibujar
  p <- ggplot() +
    geom_line(data = historico, aes(x = .data[[t_var]], y = .data[[y_var]])) +
    geom_ribbon(data = fc_ci,
                aes(x = .data[[t_var]], ymin = `95%_lower`, ymax = `95%_upper`),
                fill = "steelblue", alpha = 0.2) +
    geom_ribbon(data = fc_ci,
                aes(x = .data[[t_var]], ymin = `80%_lower`, ymax = `80%_upper`),
                fill = "steelblue", alpha = 0.35) +
    geom_line(data = fc_ci, aes(x = .data[[t_var]], y = .mean),
              color = "steelblue") +
    labs(y = y_var)
  if (!is.null(fit)) {
    ajuste <- augment(fit) |>
      filter(.data[[t_var]] >= min(historico[[t_var]]))
    p <- p + geom_line(data = ajuste, aes(x = .data[[t_var]], y = .fitted),
                       color = "firebrick", linetype = "dashed", linewidth = 0.6)
  }
  p
}


# PREDICCIÓN DE UNA SERIE COMBINANDO STL Y ARIMA.
# Estrategia general, válida tanto para series con estacionalidad como sin ella:
#   1. se descompone la serie con STL en tendencia + estacionalidad + remainder
#   2. se predice la TENDENCIA con un modelo ARIMA (al ser una curva suave, la
#      predicción no arrastra las fluctuaciones erráticas del corto plazo)
#   3. se replica el último ciclo estacional observado sobre el horizonte
#   4. el intervalo de confianza combina las DOS fuentes de incertidumbre, que
#      son independientes y por tanto suman varianzas:
#         - el error proyectado por el ARIMA, que crece con el horizonte
#         - la dispersión de la serie alrededor de su descomposición (remainder),
#           aproximadamente constante
#      de modo que sigma_total = sqrt(sigma_ARIMA^2 + sigma_remainder^2).
# Si la serie no admite estacionalidad (por ejemplo, una serie anual), STL
# devuelve solo tendencia y remainder, y el paso 3 no añade nada: el mismo código
# sirve para los dos casos.
#
#   ts                     : tsibble con UNA sola serie (sin variable key)
#   h                      : horizonte de predicción
#   periodos_atras         : nº de observaciones más recientes usadas para
#                            ajustar el modelo (Inf, por omisión, usa toda la
#                            serie). Sirve para descartar un pasado lejano que ya
#                            no representa el comportamiento actual de la serie:
#                            un tramo inicial plano o con un régimen distinto
#                            puede arrastrar la tendencia ajustada, y recortarlo
#                            deja que el modelo se centre en los años recientes.
#   ventana_tendencia      : ventana de suavizado de la tendencia (NULL = la que
#                            elige STL por defecto). A mayor valor, más lisa.
#   ventana_estacionalidad : ventana de la estacionalidad, en ciclos; "periodic"
#                            fuerza un patrón estacional constante.
# Devuelve el tsibble ampliado con las filas de la predicción y las columnas
# prediccion, lower y upper. En el tramo observado, lower y upper coinciden con
# la predicción (banda de anchura nula), porque allí no hay incertidumbre futura.
# Descomposición STL de una serie con las dos ventanas de regularidad. Función
# interna, compartida por aedv_prediccion_stl_arima() y aedv_tendencia_stl() para
# que las dos descompongan exactamente igual.
.aedv_componentes_stl <- function(ts, ventana_tendencia = NULL,
                                  ventana_estacionalidad = 13) {
  val <- tsibble::measured_vars(ts)[1]
  m   <- tsibble::guess_frequency(ts[[tsibble::index_var(ts)]])

  # SERIES SIN ESTACIONALIDAD POSIBLE (periodo 1: una serie anual).
  # La implementación de STL de feasts IGNORA en silencio la ventana de la
  # tendencia cuando no hay componente estacional: medido sobre el PIB de España,
  # las ventanas 5, 7, 13, 17, 25 y 51 devuelven exactamente la misma tendencia.
  # No es un fallo de este fichero, pero deja sin efecto un argumento que el
  # usuario cree estar usando, que es la peor clase de error: silencioso.
  #
  # Así que aquí la tendencia se calcula a mano con el MISMO suavizador que STL
  # emplea por dentro para su componente de tendencia -una regresión local lineal,
  # loess-, con el span equivalente a la ventana pedida: una ventana de w
  # observaciones sobre una serie de n es un span de w/n.
  #
  # Dos detalles medidos, no supuestos:
  #  - Con menos de 4 observaciones en la ventana loess no puede ajustar la recta
  #    local y falla. Por debajo de ese límite la única lectura sensata de la
  #    ventana es "no suavices", y eso es exactamente lo que significa la opción
  #    "sin suavizado": la tendencia ES la serie, y el resto queda a cero.
  #  - Se usa la familia gaussiana y no la robusta ("symmetric", que sería el
  #    equivalente del robust = TRUE de las llamadas a STL de abajo) porque la
  #    robusta necesita 6 observaciones por ventana y obligaría a agrandar en
  #    silencio las ventanas pequeñas, que es justo el pecado que se viene a
  #    corregir: más vale un suavizado no robusto que hace lo que se le pide que
  #    uno robusto que ignora el parámetro.
  if (m == 1 && !is.null(ventana_tendencia)) {
    idx     <- tsibble::index_var(ts)
    x       <- as.numeric(ts[[idx]])
    y       <- ts[[val]]
    ventana <- as.numeric(ventana_tendencia)

    tendencia <- if (ventana < 4) {
      y
    } else {
      as.numeric(stats::fitted(stats::loess(
        y ~ x, span = min(1, ventana / length(y)), degree = 1,
        na.action = stats::na.exclude)))
    }

    # Se devuelve con los mismos nombres de columna que components(), y SIN
    # ninguna columna "season...": es lo que mira aedv_prediccion_stl_arima()
    # para saber que no hay estacionalidad que replicar en el horizonte.
    return(tibble::tibble(
      "{idx}" := ts[[idx]],
      "{val}" := y,
      trend     = tendencia,
      remainder = y - tendencia
    ))
  }

  # Fórmula de STL: se añade season() solo si la serie admite estacionalidad
  partes <- character(0)
  if (!is.null(ventana_tendencia)) {
    partes <- c(partes, paste0("trend(window = ", ventana_tendencia, ")"))
  }
  if (m > 1) {
    ventana <- if (is.character(ventana_estacionalidad)) {
      paste0('"', ventana_estacionalidad, '"')
    } else ventana_estacionalidad
    partes <- c(partes, paste0("season(window = ", ventana, ")"))
  }
  # Si no se pide ninguna ventana y la serie no admite estacionalidad (serie
  # anual con ventana_tendencia = NULL) no hay fórmula que construir: hay que
  # llamar a STL con la variable a secas, porque una fórmula del tipo "y ~ 1"
  # STL la interpreta como un regresor externo y da error.
  if (length(partes) > 0) {
    formula_stl <- paste0(val, " ~ ", paste(partes, collapse = " + "))
    ts |>
      model(STL(as.formula(formula_stl), robust = TRUE)) |>
      components()
  } else {
    ts |>
      model(STL(!!rlang::sym(val), robust = TRUE)) |>
      components()
  }
}


# Devuelve solo la TENDENCIA que estima STL, en un tibble con el índice de la
# serie y la columna "tendencia". Es la componente que de verdad cambia al mover
# las ventanas de la descomposición: la suma de tendencia y estacionalidad
# reconstruye la serie casi igual con cualquier ventana —lo que se le quita a una
# componente reaparece en otra—, así que dibujar la tendencia es lo que permite
# ver el efecto del suavizado.
aedv_tendencia_stl <- function(ts, ventana_tendencia = NULL,
                               ventana_estacionalidad = 13) {
  idx <- tsibble::index_var(ts)
  .aedv_componentes_stl(ts, ventana_tendencia, ventana_estacionalidad) |>
    as_tibble() |>
    select(all_of(idx), tendencia = trend)
}


# EL MODELO QUE SE HA AJUSTADO, para poder enseñarlo.
# `ARIMA()` elige sola la forma del modelo, así que ni quien lo ejecuta sabe de
# antemano qué ha ajustado. Esto extrae del mable su nombre en notación
# ARIMA(p,d,q)(P,D,Q)[m] y la tabla de coeficientes estimados.
#
# Las funciones de predicción de abajo lo dejan colgado del resultado como
# ATRIBUTO y no como columna: una columna más cambiaría lo que devuelven, que
# está documentado en el capítulo 4 del libro ("tres columnas nuevas:
# prediccion, lower y upper"). El precio es que los verbos de dplyr descartan
# los atributos, así que hay que leerlo del objeto TAL COMO sale de la función,
# antes de filtrarlo o recortarlo.
.aedv_info_modelo <- function(fit) {
  list(
    descripcion = tryCatch(format(fit[[1]][[1]]), error = function(e) NA_character_),
    coeficientes = tryCatch(
      fabletools::tidy(fit) |>
        transmute(Parámetro     = term,
                  Estimación    = estimate,
                  `Error típico` = std.error,
                  `p-valor`      = p.value),
      error = function(e) NULL)
  )
}

# Accesores de lo que las funciones de predicción cuelgan del resultado
aedv_modelo <- function(resultado) attr(resultado, "aedv_modelo")

# La tendencia de STL, observada Y predicha: un tibble con el índice, la columna
# "tendencia" y la marca lógica "predicha". Solo lo trae el resultado de
# aedv_prediccion_stl_arima(), que es donde se ajusta el ARIMA de la tendencia.
aedv_tendencia <- function(resultado) attr(resultado, "aedv_tendencia")

# Desviación típica de los residuos del ajuste (en las unidades de la serie):
# sd(observado - ajuste) sobre el tramo con el que se ha entrenado. Sirve para
# comparar el ajuste ENTRE métodos de clases distintas, donde el AICc no vale
# (sección "Comparación de modelos" del capítulo 4 del libro). En
# aedv_prediccion_arima() son los residuos del modelo (columna .resid de
# augment()). En aedv_prediccion_stl_arima() NO basta con la dispersión del
# remainder: la tendencia observada la calcula STL viendo TODA la serie (pasado
# y futuro), así que "observado - trend - estacionalidad" no es un error de
# predicción real, sino que subestimaría la incertidumbre del método. Por eso
# se combinan, igual que en el intervalo de confianza de la predicción (ver la
# sección "Predicción combinando STL y ARIMA"), las dos fuentes de error
# independientes: el residuo del ARIMA que sí predice la tendencia paso a paso,
# y la dispersión del remainder, sumando varianzas.
aedv_sigma <- function(resultado) attr(resultado, "aedv_sigma")


# Escala del MASE: MAE del método trivial estacional (y_t - y_{t-m}) sobre la
# serie y. m = 1 si la serie no tiene estacionalidad (serie anual). Compartida
# por aedv_arima_accuracy() y la validación cruzada, para que las dos escalen
# con la misma regla.
.aedv_escala_mase <- function(y, m) {
  mean(abs(diff(y, lag = m)), na.rm = TRUE)
}

# MÉTRICAS DEL AJUSTE (dentro de muestra) de aedv_prediccion_arima(): RMSE, MAE
# y MASE de los residuos sobre el tramo con el que se ha entrenado el modelo.
# Dan lo mismo que accuracy() de fable sobre el modelo ajustado, pero sin tener
# que volver a ajustarlo. Miden lo bien que el modelo reproduce la serie que ya
# ha visto, NO cómo predice datos nuevos (para eso está aedv_cv_arima()).
#
# Se calculan con las columnas del propio resultado: en el tramo observado,
# "prediccion" es el ajuste del modelo, así que el residuo es observado -
# prediccion. Las filas sin ajuste (las que dejó fuera periodos_atras) y el
# tramo futuro no cuentan, y el MASE se escala sobre ese mismo tramo de
# entrenamiento.
#
# Solo vale para ARIMA: en aedv_prediccion_stl_arima() observado - prediccion
# es el remainder de una descomposición que ha visto toda la serie, y daría un
# ajuste engañosamente bueno (ver aedv_sigma()).
aedv_arima_accuracy <- function(resultado) {
  if (!is.null(attr(resultado, "aedv_tendencia"))) {
    stop("aedv_arima_accuracy() solo admite el resultado de aedv_prediccion_arima().")
  }
  idx <- tsibble::index_var(resultado)
  val <- setdiff(names(resultado), c(idx, "prediccion", "lower", "upper"))[1]

  ajuste <- resultado |>
    as_tibble() |>
    filter(!is.na(.data[[val]]), !is.na(prediccion))
  e <- ajuste[[val]] - ajuste$prediccion
  m <- max(1, round(tsibble::guess_frequency(resultado[[idx]])))

  tibble(
    RMSE = sqrt(mean(e^2)),
    MAE  = mean(abs(e)),
    MASE = mean(abs(e)) / .aedv_escala_mase(ajuste[[val]], m)
  )
}


aedv_prediccion_stl_arima <- function(ts, h, periodos_atras = Inf,
                                      ventana_tendencia = NULL,
                                      ventana_estacionalidad = 13) {
  idx <- tsibble::index_var(ts)
  val <- tsibble::measured_vars(ts)[1]
  m   <- tsibble::guess_frequency(ts[[idx]])   # 1 = sin estacionalidad posible

  # periodos_atras filtra SOLO la serie con la que se descompone y se ajusta el
  # modelo (ts_ajuste); el tramo observado que se devuelve conserva siempre la
  # serie ORIGINAL completa (ts) — ver la nota de aedv_prediccion_arima().
  ts_ajuste <- if (is.finite(periodos_atras)) utils::tail(ts, periodos_atras) else ts

  componentes <- .aedv_componentes_stl(ts_ajuste, ventana_tendencia, ventana_estacionalidad)

  # 1) tendencia y dispersión del remainder
  sd_remainder <- sd(componentes$remainder, na.rm = TRUE)

  # 2) ARIMA sobre la tendencia
  serie_tendencia <- componentes |>
    as_tibble() |>
    select(all_of(idx), trend) |>
    as_tsibble(index = all_of(idx))
  fit_tendencia <- serie_tendencia |> model(ARIMA(trend))
  fc <- fit_tendencia |> forecast(h = h)
  sd_arima <- fc |>
    hilo(level = 95) |>
    unpack_hilo(`95%`) |>
    as_tibble() |>
    mutate(.sd = (`95%_upper` - .mean) / qnorm(0.975)) |>
    pull(.sd)

  # Residuo PROPIO del ARIMA al predecir la tendencia paso a paso (no el error
  # de predicción a horizonte h de arriba, que crece con h): es la pieza que
  # falta para que aedv_sigma() sea comparable a la de ETS/ARIMA (ver más abajo)
  sigma_arima_tendencia <- fit_tendencia |> augment() |> as_tibble() |>
    pull(.resid) |> sd(na.rm = TRUE)

  # 3) réplica del último ciclo estacional completo (vacía si no hay estacionalidad)
  col_season <- setdiff(grep("^season", names(componentes), value = TRUE),
                        "season_adjust")
  estacionalidad <- if (length(col_season) == 1) {
    rep_len(componentes |> as_tibble() |> tail(m) |> pull(col_season), h)
  } else rep(0, h)

  # 4) intervalo de confianza combinando las dos varianzas
  sd_total <- sqrt(sd_arima^2 + sd_remainder^2)

  # ajuste del método, solo en el tramo con el que se ha entrenado (ts_ajuste)
  ajuste <- componentes |>
    as_tibble() |>
    transmute(
      "{idx}" := .data[[idx]],
      prediccion = trend + (if (length(col_season) == 1) .data[[col_season]] else 0),
      lower = prediccion,
      upper = prediccion
    )

  # tramo observado: la serie ORIGINAL completa, con el ajuste anterior
  # incorporado donde exista (NA en el resto, si periodos_atras recortó el
  # entrenamiento) — así la serie devuelta nunca depende de periodos_atras
  historico <- ts |>
    as_tibble() |>
    transmute("{idx}" := .data[[idx]], "{val}" := .data[[val]]) |>
    left_join(ajuste, by = idx)

  futuro <- fc |>
    as_tibble() |>
    transmute(
      "{idx}" := .data[[idx]],
      prediccion = .mean + estacionalidad,
      lower = prediccion - qnorm(0.975) * sd_total,
      upper = prediccion + qnorm(0.975) * sd_total
    )

  resultado <- bind_rows(historico, futuro) |>
    as_tsibble(index = all_of(idx))

  # La TENDENCIA, observada y predicha. La parte predicha es la que de verdad
  # extrapola este método: el ARIMA se ajusta sobre la tendencia, no sobre la
  # serie, y la predicción final es esta curva más el ciclo estacional replicado.
  # Dibujarla es lo que permite ver de dónde sale la predicción.
  attr(resultado, "aedv_tendencia") <- bind_rows(
    componentes |> as_tibble() |>
      transmute("{idx}" := .data[[idx]], tendencia = trend, predicha = FALSE),
    fc |> as_tibble() |>
      transmute("{idx}" := .data[[idx]], tendencia = .mean, predicha = TRUE)
  )
  attr(resultado, "aedv_modelo") <- .aedv_info_modelo(fit_tendencia)
  # Dos fuentes de error INDEPENDIENTES, igual que en el intervalo de
  # confianza de la predicción: el ajuste del ARIMA a la tendencia y la
  # dispersión del remainder. Se combinan sumando varianzas.
  attr(resultado, "aedv_sigma")  <- sqrt(sigma_arima_tendencia^2 + sd_remainder^2)
  resultado
}


# PREDICCIÓN CON UN MODELO ARIMA, devuelta en el mismo formato que
# aedv_prediccion_stl_arima() para poder comparar ambos métodos con el mismo
# código. En el tramo observado la columna "prediccion" contiene los valores
# AJUSTADOS por el modelo (fitted), es decir, lo que el modelo dice que valía la
# serie en cada instante; en el tramo futuro, la predicción y su intervalo al 95%.
# periodos_atras restringe el AJUSTE (no la serie devuelta) a las observaciones
# más recientes: Inf, por omisión, usa toda la serie. Filtra localmente la serie
# con la que se entrena el modelo (ts_ajuste); el tramo observado que se
# devuelve conserva siempre la serie ORIGINAL completa que se ha recibido, con
# el ajuste incorporado solo donde exista (NA en el resto). Los argumentos
# adicionales (...) se pasan a ARIMA(), lo que permite por ejemplo afinar la
# búsqueda con stepwise = FALSE, approx = FALSE.
aedv_prediccion_arima <- function(ts, h, periodos_atras = Inf, ...) {
  idx <- tsibble::index_var(ts)
  val <- tsibble::measured_vars(ts)[1]

  ts_ajuste <- if (is.finite(periodos_atras)) utils::tail(ts, periodos_atras) else ts

  fit <- ts_ajuste |> model(ARIMA(!!rlang::sym(val), ...))
  fc  <- fit |>
    forecast(h = h) |>
    hilo(level = 95) |>
    unpack_hilo(`95%`)

  # ajuste del modelo, solo en el tramo con el que se ha entrenado (ts_ajuste)
  augmentado <- fit |> augment() |> as_tibble()
  ajuste <- augmentado |>
    transmute("{idx}" := .data[[idx]], prediccion = .fitted, lower = .fitted, upper = .fitted)

  # tramo observado: la serie ORIGINAL completa, con el ajuste anterior
  # incorporado donde exista
  historico <- ts |>
    as_tibble() |>
    transmute("{idx}" := .data[[idx]], "{val}" := .data[[val]]) |>
    left_join(ajuste, by = idx)

  # tramo futuro: la predicción y su intervalo de confianza
  futuro <- fc |>
    as_tibble() |>
    transmute(
      "{idx}" := .data[[idx]],
      prediccion = .mean,
      lower = `95%_lower`,
      upper = `95%_upper`
    )

  resultado <- bind_rows(historico, futuro) |>
    as_tsibble(index = all_of(idx))
  attr(resultado, "aedv_modelo") <- .aedv_info_modelo(fit)
  attr(resultado, "aedv_sigma")  <- sd(augmentado$.resid, na.rm = TRUE)
  resultado
}


# Dibuja el resultado de aedv_prediccion_arima(): la serie observada en gris, el
# AJUSTE del modelo sobre esa misma serie en rojo discontinuo, y la predicción en
# azul con su banda del intervalo de confianza. Ver el ajuste superpuesto es lo
# que permite juzgar si el modelo reproduce bien la serie con la que se entrenó.
# Colores y etiquetas compartidos por aedv_dibujar_arima() y
# aedv_dibujar_stl_arima(): mapear "colour" a estas etiquetas (en vez de fijar
# el color fuera de aes()) es lo que permite que aedv_dibujar_..._dinamico()
# etiquete cada traza correctamente en la leyenda y en el cuadro emergente de
# plotly, sin que la versión estática (con la leyenda oculta) cambie de aspecto.
.AEDV_COLORES_PREDICCION <- c(
  "Observado"              = "grey45",
  "Ajuste del modelo"      = "firebrick",
  "Predicción"             = "steelblue",
  "Límite inferior (95%)"  = "steelblue",
  "Límite superior (95%)"  = "steelblue"
)

# Capas comunes a aedv_dibujar_arima() y aedv_dibujar_stl_arima(): la banda
# sombreada, la serie observada, lo que el método dice del tramo conocido
# ("Observado"/"Ajuste del modelo" en el ARIMA directo; "Descomposición" en
# STL + ARIMA, mismo color) y la predicción con sus dos límites explícitos.
.aedv_capas_prediccion <- function(resultado, idx, val, etiqueta_ajuste) {
  observado <- resultado |> as_tibble() |> filter(!is.na(.data[[val]]))
  predicho  <- resultado |> as_tibble() |> filter(is.na(.data[[val]]))
  colores   <- .AEDV_COLORES_PREDICCION
  names(colores)[2] <- etiqueta_ajuste

  list(
    geom_ribbon(data = resultado, aes(x = .data[[idx]], ymin = lower, ymax = upper),
                fill = "steelblue", alpha = 0.2),
    # la serie observada va debajo y con más grosor, para que no la tape el ajuste
    geom_line(data = observado,
              aes(x = .data[[idx]], y = .data[[val]], colour = "Observado"),
              linewidth = 0.9),
    geom_line(data = observado,
              aes(x = .data[[idx]], y = prediccion, colour = !!etiqueta_ajuste),
              linetype = "dashed", linewidth = 0.45),
    geom_line(data = predicho,
              aes(x = .data[[idx]], y = prediccion, colour = "Predicción"),
              linewidth = 1),
    # límites del intervalo, finos y punteados: solo se ven bien diferenciados
    # de la banda sombreada al pasar el ratón, que es donde de verdad importan
    geom_line(data = predicho,
              aes(x = .data[[idx]], y = lower, colour = "Límite inferior (95%)"),
              linewidth = 0.6, linetype = "dotted"),
    geom_line(data = predicho,
              aes(x = .data[[idx]], y = upper, colour = "Límite superior (95%)"),
              linewidth = 0.6, linetype = "dotted"),
    scale_colour_manual(values = colores, breaks = names(colores)),
    theme(legend.position = "none")
  )
}


aedv_dibujar_arima <- function(resultado,
                               titulo    = NULL,
                               subtitulo = "Ajuste del modelo en rojo discontinuo; predicción y banda al 95% en azul",
                               x_lab     = "",
                               y_lab     = NULL) {
  idx <- tsibble::index_var(resultado)
  val <- setdiff(names(resultado), c(idx, "prediccion", "lower", "upper"))[1]
  if (is.null(y_lab)) y_lab <- val

  ggplot() +
    .aedv_capas_prediccion(resultado, idx, val, "Ajuste del modelo") +
    labs(x = x_lab, y = y_lab, title = titulo, subtitle = subtitulo, colour = NULL)
}


# ¿Qué TENDENCIA de STL se dibuja cuando no se pide ninguna en concreto? La que
# viaja como atributo del resultado (ver aedv_tendencia()). La usan por igual
# aedv_dibujar_stl_arima() y su versión dinámica, para que las dos enseñen lo
# mismo sin que haya que pasarles nada.
#
# Devuelve NULL -no se dibuja tendencia- en dos casos:
#   - Si el resultado ha perdido el atributo "aedv_tendencia", por ejemplo al
#     filtrarlo con un verbo de dplyr (ver la nota sobre atributos al principio
#     de aedv_prediccion_stl_arima()).
#   - Si la serie NO tiene estacionalidad (una serie anual, como el PIB de un
#     país): en ese caso la predicción del método ES la tendencia extrapolada
#     -la estacionalidad que se le suma vale cero-, así que la línea naranja
#     quedaría exactamente debajo de la roja y de la azul, sin enseñar nada
#     nuevo.
#
# Un `tendencia` pasado a mano NO pasa por aquí: se dibuja tal cual llegue. Es
# lo que necesita un cuadro de mando que la recalcule al mover las ventanas del
# suavizado (ver Dashboard6.Rmd), donde la línea naranja es justo lo que hace
# visible el efecto de los controles aunque la serie sea anual.
.aedv_tendencia_por_omision <- function(resultado) {
  tendencia <- attr(resultado, "aedv_tendencia")
  if (is.null(tendencia)) return(NULL)
  idx <- tsibble::index_var(resultado)
  if (tsibble::guess_frequency(resultado[[idx]]) == 1) return(NULL)
  tendencia
}


# Capa de la TENDENCIA de STL, observada y predicha (ver aedv_tendencia()): la
# curva SIN estacionalidad de la que sale la predicción de aedv_dibujar_stl_arima(),
# distinta de "Descomposición" (que sí lleva sumada la estacionalidad) y de
# "Predicción" (que además incluye el intervalo de confianza). Se dibuja en
# naranja, sólida en el tramo observado y discontinua en el predicho,
# arrastrando el último punto observado para que las dos líneas queden
# empalmadas sin un salto visual (mismo criterio que usa, por separado, la
# versión dinámica de este gráfico).
.aedv_capa_tendencia <- function(tendencia, idx) {
  if (is.null(tendencia)) return(list())

  tendencia <- as_tibble(tendencia)
  observada <- tendencia |> filter(!predicha)
  prevista  <- tendencia |> filter(predicha)
  if (nrow(prevista) > 0 && nrow(observada) > 0) {
    prevista <- bind_rows(tail(observada, 1), prevista)
  }

  list(
    geom_line(data = observada, aes(x = .data[[idx]], y = tendencia),
              colour = "darkorange", linewidth = 0.6),
    geom_line(data = prevista, aes(x = .data[[idx]], y = tendencia),
              colour = "darkorange", linewidth = 0.6, linetype = "dashed")
  )
}


# Dibuja el resultado de aedv_prediccion_stl_arima(): la serie observada en gris,
# la tendencia en naranja (cuando la serie tiene estacionalidad, ver
# .aedv_tendencia_por_omision()), la predicción en azul y la banda del intervalo
# de confianza sombreada.
#
# El argumento "tendencia" admite el tibble que devuelve aedv_tendencia_stl(),
# con las columnas del índice, `tendencia` y `predicha`; por omisión se toma la
# que trae el propio resultado. Con `tendencia = NULL` la línea naranja no se
# dibuja: es lo que hace aedv_dibujar_stl_arima_dinamico(), que la añade después
# como traza de plotly.
aedv_dibujar_stl_arima <- function(resultado,
                                        titulo    = NULL,
                                        subtitulo = NULL,
                                        x_lab     = "",
                                        y_lab     = NULL,
                                        tendencia = .aedv_tendencia_por_omision(resultado)) {
  idx <- tsibble::index_var(resultado)
  val <- setdiff(names(resultado), c(idx, "prediccion", "lower", "upper"))[1]
  if (is.null(y_lab)) y_lab <- val

  capa_tendencia <- .aedv_capa_tendencia(tendencia, idx)
  if (is.null(subtitulo)) {
    subtitulo <- if (length(capa_tendencia) == 0)
      "Descomposición en rojo discontinuo; predicción y banda al 95% en azul"
    else
      "Tendencia en naranja; descomposición en rojo discontinuo; predicción y banda al 95% en azul"
  }
  # el texto de este subtítulo (sobre todo con la tendencia) es más largo que el
  # de aedv_dibujar_arima() y ggplot no lo envuelve solo: partirlo en varias
  # líneas evita que se salga del ancho de la imagen
  subtitulo <- stringr::str_wrap(subtitulo, width = 50)

  ggplot() +
    # mismo código de colores que aedv_dibujar_arima(): la serie observada en
    # gris, lo que el modelo dice del pasado en rojo discontinuo y la predicción
    # en azul
    .aedv_capas_prediccion(resultado, idx, val, "Descomposición") +
    capa_tendencia +
    labs(x = x_lab, y = y_lab, title = titulo, subtitle = subtitulo, colour = NULL)
}


# VALIDACIÓN CRUZADA TEMPORAL con ventana deslizante para las funciones de
# predicción de este fichero. Sirve para lo que el ajuste sobre los datos
# conocidos no puede decir: cómo se comporta el método sobre datos que NO ha
# visto.
#
# El procedimiento es el de la sección de validación cruzada del capítulo de
# series temporales: se toma una ventana de entrenamiento (el PLIEGUE), se
# predice h periodos, se comparan las predicciones con lo realmente observado, y
# la ventana avanza UNA observación en el tiempo. Nunca se entrena con datos
# posteriores al punto que se predice. El número de pliegues no se elige: se
# usan TODOS los orígenes de predicción que caben en la serie —el máximo detalle
# posible—, así que se deduce de `ts`, `tamano_pliegue` y `h`.
#
# Hay una función por método, `aedv_cv_arima()` y `aedv_cv_stl_arima()`, en vez
# de una única función genérica con la función de predicción como parámetro: las
# dos comparten esta implementación privada.
#
#   ts             : tsibble con UNA sola serie
#   h              : horizonte de predicción de cada pliegue
#   predictor      : función de predicción interna, con la firma predictor(ts, h)
#   tamano_pliegue : nº de observaciones con las que se entrena en cada pliegue.
#                    Cada pliegue reajusta el modelo, así que cuanto más pequeño
#                    es (para una serie dada), más pliegues salen y más tarda.
#   nivel          : nivel del intervalo de confianza
#   banda          : cómo se estima la anchura de la banda a cada horizonte.
#                    "tipica" (por omisión): la predicción más/menos el error
#                    típico observado a ese horizonte (su RMSE), escalado por el
#                    cuantil t que corresponde al nº de pliegues. "empirica": los
#                    cuantiles del error, sin suponer simetría, como en el
#                    capítulo de series temporales del libro; solo tiene sentido
#                    con MUCHOS pliegues, porque con pocos errores el cuantil al
#                    2,5 % es poco más que el mínimo observado y la banda se
#                    ensancha a saltos.
#
# Devuelve una lista con:
#   metricas   : RMSE, MAE, MAPE, MASE y MdAE fuera de muestra, sobre todos los
#                pliegues. El MASE se escala con el MAE del método trivial
#                estacional (y_t - y_{t-m}) calculado sobre toda la serie ts
#   bandas     : por cada horizonte hh, el RMSE, MAE y MASE de ese horizonte, con
#                cuántos errores se han estimado, y los desplazamientos que hay
#                que sumar a la predicción para obtener la banda. No dependen de
#                los supuestos del modelo, cuyos propios intervalos suelen quedar
#                demasiado estrechos.
#   n_pliegues, tamano_pliegue : la configuración realmente aplicada (n_pliegues
#                deducido, no elegido)
.aedv_cv_temporal <- function(ts, h, predictor, tamano_pliegue,
                              nivel = 0.95, banda = c("tipica", "empirica")) {
  banda <- match.arg(banda)
  idx <- tsibble::index_var(ts)
  val <- tsibble::measured_vars(ts)[1]
  n   <- nrow(ts)
  periodo <- tsibble::guess_frequency(ts[[idx]])

  # Escala del MASE: MAE del método trivial estacional (y_t - y_{t-m}), sobre la
  # serie completa. m = 1 si no hay estacionalidad (periodo anual o guess_frequency
  # devuelve 1), igual que en la sección de validación cruzada del libro.
  m <- max(1, round(periodo))
  escala_mase <- .aedv_escala_mase(ts[[val]], m)

  if (tamano_pliegue + h > n) {
    stop("No hay sitio para un solo pliegue: el tamaño del pliegue (",
         tamano_pliegue, ") más el horizonte (", h, ") superan las ",
         n, " observaciones de la serie.")
  }

  # Todos los orígenes de predicción posibles, avanzando de uno en uno: el
  # número de pliegues es consecuencia del tamaño de la serie, del pliegue y del
  # horizonte, no una elección aparte.
  origenes <- seq(tamano_pliegue, n - h)

  errores <- vector("list", length(origenes))
  for (k in seq_along(origenes)) {
    o <- origenes[k]
    entrenamiento <- ts[(o - tamano_pliegue + 1):o, ]

    prediccion <- predictor(entrenamiento, h) |>
      as_tibble() |>
      filter(is.na(.data[[val]])) |>
      mutate(hh = dplyr::row_number()) |>
      select(hh, prediccion)

    # Valores realmente observados después del origen, con su horizonte
    real <- ts[(o + 1):min(o + h, n), ] |>
      as_tibble() |>
      mutate(hh = dplyr::row_number()) |>
      select(hh, observado = all_of(val))

    errores[[k]] <- inner_join(real, prediccion, by = "hh") |>
      mutate(pliegue = k, error = observado - prediccion)
  }
  errores <- bind_rows(errores)

  alfa <- (1 - nivel) / 2
  bandas <- errores |>
    summarise(
      n_errores = dplyr::n(),
      # error típico a este horizonte: incluye tanto el sesgo como la dispersión
      rmse = sqrt(mean(error^2)),
      mae  = mean(abs(error)),
      mase = mean(abs(error)) / escala_mase, # mismo denominador que el MASE global
      q_inferior_emp = quantile(error, alfa),
      q_superior_emp = quantile(error, 1 - alfa),
      .by = hh
    ) |>
    arrange(hh)

  if (banda == "tipica") {
    # Se usa el cuantil de la t de Student, no el de la normal, porque el error
    # típico está estimado con muy pocos pliegues: con 6 errores, t = 2,57 frente
    # al 1,96 de la normal. Con un solo error no hay dispersión que estimar.
    bandas <- bandas |>
      mutate(
        t = ifelse(n_errores > 1, qt(1 - alfa, df = n_errores - 1), NA_real_),
        q_inferior = -t * rmse,
        q_superior =  t * rmse
      ) |>
      select(hh, n_errores, rmse, mae, mase, q_inferior, q_superior)
  } else {
    bandas <- bandas |>
      mutate(q_inferior = q_inferior_emp, q_superior = q_superior_emp) |>
      select(hh, n_errores, rmse, mae, mase, q_inferior, q_superior)
  }

  list(
    metricas = tibble(
      RMSE       = sqrt(mean(errores$error^2)),
      MAE        = mean(abs(errores$error)),
      `MAPE (%)` = 100 * mean(abs(errores$error / errores$observado)),
      MASE       = mean(abs(errores$error)) / escala_mase,
      MdAE       = median(abs(errores$error))
    ),
    bandas = bandas,
    n_pliegues     = length(origenes),
    tamano_pliegue = tamano_pliegue
  )
}

# Validación cruzada temporal de aedv_prediccion_arima(). Por omisión banda =
# "empirica": es habitual querer después los cuantiles empíricos por horizonte
# para aedv_aplicar_bandas_cv() (ver la sección de intervalos no paramétricos
# del capítulo), y así el mismo resultado sirve para las métricas y para las
# bandas sin recalcular nada. Los argumentos adicionales (...) se pasan a
# ARIMA(), igual que en aedv_prediccion_arima().
aedv_cv_arima <- function(ts, h, tamano_pliegue, nivel = 0.95,
                          banda = "empirica", ...) {
  .aedv_cv_temporal(ts, h,
                    predictor      = function(x, h) aedv_prediccion_arima(x, h, ...),
                    tamano_pliegue = tamano_pliegue, nivel = nivel, banda = banda)
}

# Validación cruzada temporal de aedv_prediccion_stl_arima(). banda = "empirica"
# por omisión, igual que aedv_cv_arima() y por el mismo motivo.
aedv_cv_stl_arima <- function(ts, h, tamano_pliegue, nivel = 0.95,
                              banda = "empirica",
                              ventana_tendencia = NULL,
                              ventana_estacionalidad = 13) {
  .aedv_cv_temporal(ts, h,
                    predictor      = function(x, h) aedv_prediccion_stl_arima(
                      x, h, ventana_tendencia = ventana_tendencia,
                      ventana_estacionalidad = ventana_estacionalidad),
                    tamano_pliegue = tamano_pliegue, nivel = nivel, banda = banda)
}


# Sustituye los intervalos de confianza del modelo por los que estiman
# aedv_cv_arima()/aedv_cv_stl_arima(): en el tramo predicho, la banda pasa a ser
# la predicción más los desplazamientos calculados con los errores observados a
# ese mismo horizonte.
aedv_aplicar_bandas_cv <- function(resultado, bandas) {
  idx <- tsibble::index_var(resultado)
  val <- setdiff(names(resultado), c(idx, "prediccion", "lower", "upper"))[1]

  observado <- resultado |> as_tibble() |> filter(!is.na(.data[[val]]))
  futuro    <- resultado |> as_tibble() |> filter(is.na(.data[[val]])) |>
    mutate(hh = dplyr::row_number()) |>
    left_join(bandas, by = "hh") |>
    mutate(
      # si el horizonte va más allá de lo evaluado por la validación cruzada, se
      # conserva la banda del modelo en vez de inventar nada
      lower = ifelse(is.na(q_inferior), lower, prediccion + q_inferior),
      upper = ifelse(is.na(q_superior), upper, prediccion + q_superior)
    ) |>
    select(-hh, -q_inferior, -q_superior)

  # bind_rows() construye un objeto nuevo y los verbos de dplyr descartan los
  # atributos, así que hay que volver a colgarlos todos: si no, al marcar la
  # validación cruzada el modelo y la tendencia desaparecerían de la pantalla
  # sin más, y aedv_sigma() dejaría de funcionar.
  nuevo <- bind_rows(observado, futuro) |> as_tsibble(index = all_of(idx))
  attr(nuevo, "aedv_modelo")    <- attr(resultado, "aedv_modelo")
  attr(nuevo, "aedv_tendencia") <- attr(resultado, "aedv_tendencia")
  attr(nuevo, "aedv_sigma")     <- attr(resultado, "aedv_sigma")
  nuevo
}


# Convierte con ggplotly() un gráfico de predicción de los dos anteriores y
# arregla lo que la conversión estropea. Función interna: las que se usan son
# aedv_dibujar_arima_dinamico() y aedv_dibujar_stl_arima_dinamico().
#
# El cuadro emergente NO lo genera ggplotly() a partir de las 5 líneas: con
# "x unified"/"compare data on hover", plotly busca en CADA traza su punto
# real más cercano al cursor sin límite de distancia (ni "hoverdistance" lo
# acota: probado a mano bajándolo hasta 1 px y el fallo seguía). Como el
# tramo histórico y la predicción son disjuntos, esa búsqueda encuentra
# igualmente el primer punto del tramo contrario y lo muestra como si
# coincidiera con la fecha señalada (p. ej. "Predicción" asomando en pleno
# histórico, con el valor de dentro de tres meses). Por eso las 5 líneas se
# dejan solo para VERSE (con su color y su entrada en la leyenda) y se les
# quita el cuadro emergente propio; dos trazas invisibles (una por tramo, ver
# más abajo) son las únicas que responden al ratón, cada una con su propio
# texto ya armado en R con los campos que de verdad le corresponden a esa
# fecha, sin que plotly tenga que adivinar nada.
#
# FORMATO DE LOS NÚMEROS DEL CUADRO EMERGENTE, elegido según la magnitud de la
# serie. Un formato fijo sin decimales va bien para una serie de cientos de miles
# de pasajeros, pero convierte un PIB de 10,74 en un escueto "11", y con él las
# tres líneas del cuadro emergente acaban mostrando el mismo número.
#
# Los decimales se deciden por la magnitud de los valores PEQUEÑOS, no por la del
# máximo: una serie larga puede recorrer varios órdenes de magnitud —el PIB de
# España va de 10,7 en 1960 a 1.631,9 en 2023—, y quien mire el arranque de la
# serie es justo el que necesita los decimales. Al mayor le sobran, pero un
# decimal de más no engaña a nadie; uno de menos, sí.
#
# Se toma el PERCENTIL 10 y no el mínimo para que un único valor cercano a cero
# (el cruce de una serie que cambia de signo, o un extremo de la banda del
# intervalo) no llene de decimales toda la serie. Y se limita a 4 decimales.
.aedv_formato_numero <- function(valores) {
  valores  <- abs(valores[is.finite(valores) & valores != 0])
  if (length(valores) == 0) return(",.2f")
  magnitud <- stats::quantile(valores, 0.1, names = FALSE)
  paste0(",.", max(0, min(4, 3 - floor(log10(magnitud)))), "f")
}


# El argumento "tendencia" (opcional, solo lo usa la versión de STL + ARIMA) es
# la componente de tendencia de la descomposición, que se superpone en naranja.
# Se añade AQUÍ, y no fuera con un add_lines() sobre el resultado, porque la "x"
# de esta función ya no es la que trae el gráfico: es texto ISO. Una traza
# añadida por fuera con el índice tal cual (un yearmonth, por ejemplo) se
# serializa como NÚMERO de días, y bajo un eje "type = date" plotly lee ese
# número como MILISEGUNDOS: la línea se va al 1 de enero de 1970 y, con
# autorange, arrastra el eje entero hasta allí dejando la serie aplastada contra
# el borde derecho.
.aedv_prediccion_a_plotly <- function(p, titulo, es_fecha, resultado, idx, val,
                                      etiqueta_ajuste, tendencia = NULL) {
  g <- plotly::ggplotly(p, tooltip = "y")

  # MISMA TRAMPA QUE LA DE ARRIBA, PERO EN LAS MARCAS DEL EJE. ggplotly no deja
  # que plotly elija las marcas: las fija a mano con tickmode = "array" y una
  # lista de posiciones calculada sobre el eje ORIGINAL, es decir, en días desde
  # 1970 (10957, 12784...). Al pasar el eje a "type = date" esos números se leen
  # como milisegundos, así que las cinco marcas se amontonan en 1970, fuera de la
  # vista: el eje horizontal se queda SIN NINGUNA etiqueta y sin sus líneas de
  # rejilla, y encima sin ningún error. Se borran para que plotly calcule las
  # marcas él mismo, que sobre un eje de fechas es justo lo que se quiere (elige
  # el paso según el zoom). En un eje numérico -una serie anual- las posiciones
  # de ggplotly son correctas y se dejan como están.
  if (es_fecha) {
    g$x$layout$xaxis$tickmode <- NULL
    g$x$layout$xaxis$tickvals <- NULL
    g$x$layout$xaxis$ticktext <- NULL
  }

  for (i in seq_along(g$x$data)) {
    relleno <- g$x$data[[i]]$fill
    if (!is.null(relleno) && relleno == "toself") g$x$data[[i]]$showlegend <- FALSE
    g$x$data[[i]]$hoverinfo <- "skip"

    # ggplot dibuja un índice yearmonth/Date convirtiéndolo primero a Date, y
    # ese Date se representa por dentro como un número (días desde 1970-01-01):
    # es lo que llega aquí, y plotly lo muestra tal cual ("18536.00") porque no
    # sabe que son días. Se convierte a fecha ISO para que lo interprete bien.
    if (es_fecha && !is.null(g$x$data[[i]]$x)) {
      g$x$data[[i]]$x <- format(as.Date(g$x$data[[i]]$x, origin = "1970-01-01"))
    }
  }

  datos <- resultado |> as_tibble()
  # misma representación de "x" que las trazas de arriba (texto ISO si hay
  # fecha de calendario), para que compartan un único eje sin ambigüedad
  fecha_x    <- if (es_fecha) format(as.Date(datos[[idx]])) else datos[[idx]]
  etiqueta_x <- if (es_fecha) format(as.Date(datos[[idx]]), "%Y-%m") else as.character(datos[[idx]])
  es_futuro  <- is.na(datos[[val]])

  historico <- tibble(x = fecha_x, etiqueta = etiqueta_x,
                       observado = datos[[val]], ajuste = datos$prediccion) |>
    filter(!es_futuro)
  futuro <- tibble(x = fecha_x, etiqueta = etiqueta_x,
                    prediccion = datos$prediccion, lower = datos$lower, upper = datos$upper) |>
    filter(es_futuro)

  # La tendencia, con la MISMA representación de "x" que el resto, y unida a los
  # dos tramos por esa "x" para que su valor pueda leerse en el cuadro emergente
  # (la línea naranja, como las demás, no lleva cuadro propio). Si trae la marca
  # "predicha" se separa en dos tramos, porque el tramo predicho de la tendencia
  # es una EXTRAPOLACIÓN del ARIMA y no debe parecer un dato observado.
  if (!is.null(tendencia)) {
    t0   <- tendencia |> as_tibble()
    tend <- tibble(
      x         = if (es_fecha) format(as.Date(t0[[idx]])) else t0[[idx]],
      tendencia = t0$tendencia,
      predicha  = if ("predicha" %in% names(t0)) t0$predicha else FALSE
    )
    historico <- historico |> left_join(tend |> select(x, tendencia), by = "x")
    futuro    <- futuro    |> left_join(tend |> select(x, tendencia), by = "x")
  }

  # Los decimales se deciden UNA VEZ, con todo lo que se va a mostrar, para que
  # las tres líneas del cuadro emergente y los dos tramos usen el mismo formato
  fmt <- .aedv_formato_numero(c(datos[[val]], datos$prediccion,
                                datos$lower, datos$upper,
                                if (!is.null(tendencia)) historico$tendencia))

  # customdata se pasa como valor literal (no como fórmula "~columna"), que es
  # la vía NSE de add_trace(): con ella, un customdata de un solo valor por
  # punto ("ajuste") llegaba bien, pero un customdata de DOS valores por punto
  # ("cbind(lower, upper)", para los límites) se perdía por completo -la traza
  # quedaba sin customdata y "%{customdata[0]}" salía literal en el rótulo-.
  # Las líneas de la tendencia van antes que las trazas invisibles, para que
  # estas sigan siendo las últimas y capten el ratón en todo el gráfico. El
  # tramo observado va continuo y el predicho discontinuo; el predicho arrastra
  # la última observación para que las dos líneas queden empalmadas y no se vea
  # un salto donde no lo hay.
  if (!is.null(tendencia)) {
    observada <- tend |> filter(!predicha)
    prevista  <- tend |> filter(predicha)
    if (nrow(prevista) > 0 && nrow(observada) > 0) {
      prevista <- bind_rows(tail(observada, 1), prevista)
    }

    g <- g |>
      plotly::add_lines(
        data = observada, x = ~x, y = ~tendencia, name = "Tendencia STL",
        line = list(color = "darkorange", width = 2),
        hoverinfo = "skip", inherit = FALSE
      )
    if (nrow(prevista) > 0) {
      g <- g |>
        plotly::add_lines(
          data = prevista, x = ~x, y = ~tendencia, name = "Tendencia predicha",
          line = list(color = "darkorange", width = 2, dash = "dash"),
          hoverinfo = "skip", inherit = FALSE
        )
    }
  }

  # AJUSTE AUSENTE: con periodos_atras el modelo se entrena solo con el tramo
  # final de la serie, y las observaciones anteriores se devuelven sin ajuste
  # (NA), igual que su tendencia. El gráfico estático lo resuelve solo -la línea
  # roja simplemente empieza más tarde-, pero aquí no basta con dejar el NA:
  # plotly formatea un customdata nulo como "0", y el cuadro emergente anunciaría
  # un ajuste de cero justo donde no hay ninguno. Se separan en dos trazas
  # invisibles: las observaciones con ajuste lo muestran, y las anteriores
  # enseñan solo el valor observado.
  sin_ajuste <- historico |> filter(is.na(ajuste))
  historico  <- historico |> filter(!is.na(ajuste))

  # Con la tendencia hay DOS valores por punto en el tramo histórico, y un
  # customdata de más de un valor por punto solo llega si se pasa como lista de
  # vectores (ver la nota de abajo sobre los límites del intervalo)
  customdata_historico <- if (is.null(tendencia)) {
    historico$ajuste
  } else {
    lapply(seq_len(nrow(historico)),
           \(i) c(historico$ajuste[i], historico$tendencia[i]))
  }
  plantilla_historico <- paste0(
    "%{text}<br>Observado: %{y:", fmt, "}<br>", etiqueta_ajuste, ": ",
    if (is.null(tendencia)) paste0("%{customdata:", fmt, "}") else
      paste0("%{customdata[0]:", fmt, "}<br>Tendencia STL: %{customdata[1]:",
             fmt, "}"),
    "<extra></extra>")

  g <- g |>
    plotly::add_trace(
      data = historico, x = ~x, y = ~observado, text = ~etiqueta,
      customdata = customdata_historico,
      type = "scatter", mode = "markers", marker = list(opacity = 0, size = 16),
      hovertemplate = plantilla_historico,
      hoverlabel = list(bgcolor = "grey45", font = list(color = "white")),
      showlegend = FALSE, inherit = FALSE
    )

  if (nrow(sin_ajuste) > 0) {
    g <- g |>
      plotly::add_trace(
        data = sin_ajuste, x = ~x, y = ~observado, text = ~etiqueta,
        type = "scatter", mode = "markers", marker = list(opacity = 0, size = 16),
        hovertemplate = paste0("%{text}<br>Observado: %{y:", fmt, "}<extra></extra>"),
        hoverlabel = list(bgcolor = "grey45", font = list(color = "white")),
        showlegend = FALSE, inherit = FALSE
      )
  }

  g <- g |>
    plotly::add_trace(
      data = futuro, x = ~x, y = ~prediccion, text = ~etiqueta,
      customdata = lapply(seq_len(nrow(futuro)), \(i) {
        if (is.null(tendencia)) c(futuro$lower[i], futuro$upper[i])
        else c(futuro$lower[i], futuro$upper[i], futuro$tendencia[i])
      }),
      type = "scatter", mode = "markers", marker = list(opacity = 0, size = 16),
      hovertemplate = paste0("%{text}<br>Predicción: %{y:", fmt, "}",
                              "<br>Límite inferior (95%): %{customdata[0]:", fmt, "}",
                              "<br>Límite superior (95%): %{customdata[1]:", fmt, "}",
                              if (!is.null(tendencia))
                                paste0("<br>Tendencia predicha: %{customdata[2]:", fmt, "}"),
                              "<extra></extra>"),
      hoverlabel = list(bgcolor = "steelblue", font = list(color = "white")),
      showlegend = FALSE, inherit = FALSE
    )

  g |>
    plotly::layout(
      # "closest" (el modo por omisión de plotly) muestra solo el punto real
      # más próximo al cursor, nunca una mezcla: es lo que hace imposible que
      # vuelva a colarse un valor del tramo contrario, y cada una de las dos
      # trazas invisibles ya trae en su plantilla todo lo que hacía falta
      # mostrar junto (el ajuste junto al observado; los dos límites junto a
      # la predicción).
      hovermode = "closest",
      # el título va aparte, sin el subtítulo: la leyenda ya explica los
      # colores, y forzarlo dentro del título lo desbordaba en gráficos
      # estrechos (el título es corto; el subtítulo, una frase entera)
      title = if (!is.null(titulo) && nzchar(titulo)) list(text = titulo),
      # HAY QUE ABRIRLE SITIO AL TÍTULO. El ggplot de partida no lleva título
      # -se pone aquí, en el layout-, así que ggplotly reserva arriba un margen
      # de andar por casa (23 px) en el que el título no cabe y sale cortado por
      # la mitad. El margen inferior se deja en manos de "automargin", que ya
      # crece solo para las etiquetas del eje; el superior no, porque el título
      # no es parte de ningún eje y automargin no lo tiene en cuenta.
      margin = if (!is.null(titulo) && nzchar(titulo)) list(t = 44),
      # el ggplot estático oculta su leyenda (theme(legend.position="none")), y
      # ggplotly copia ese "showlegend" al layout entero, por encima del de
      # cada traza: hay que reactivarlo aquí explícitamente. Se deja la
      # POSICIÓN de la leyenda en la que trae plotly por omisión (vertical, a
      # la derecha): es la única que ajusta sola el margen del gráfico para no
      # solaparse con nada; una leyenda horizontal fija arriba, con cinco
      # entradas, se envuelve a dos líneas en un panel estrecho y tapa el título.
      showlegend = TRUE,
      # autorange = TRUE es imprescindible: ggplotly calculó el rango del eje
      # con los números de días de ANTES de convertir "x" a texto ISO, y ese
      # rango viejo se interpreta como MILISEGUNDOS bajo type="date" -una
      # ventana de apenas segundos junto al 1 de enero de 1970-, dejando el
      # gráfico en blanco.
      xaxis = if (es_fecha) list(type = "date", autorange = TRUE) else list()
    )
}


# Versiones DINÁMICAS de las dos funciones de dibujo anteriores, para cuadros de
# mando y documentos HTML: mismos argumentos y mismo aspecto (salvo el
# subtítulo, que aquí no hace falta: los colores los explica la leyenda), pero
# el resultado es un gráfico de plotly con el que se puede interactuar
# (consultar valores con el ratón, ampliar un tramo arrastrando y volver a la
# vista completa con doble clic). En un documento PDF hay que usar las
# versiones estáticas.
aedv_dibujar_arima_dinamico <- function(resultado, titulo = NULL,
                                        x_lab = "", y_lab = NULL) {
  p <- aedv_dibujar_arima(resultado, x_lab = x_lab, y_lab = y_lab)
  idx <- tsibble::index_var(resultado)
  val <- setdiff(names(resultado), c(idx, "prediccion", "lower", "upper"))[1]
  .aedv_prediccion_a_plotly(p, titulo, .aedv_indice_es_fecha(resultado),
                             resultado, idx, val, "Ajuste del modelo")
}


# La TENDENCIA de la descomposición se superpone en naranja igual que en la
# versión estática, y con el mismo criterio por omisión (ver
# .aedv_tendencia_por_omision()): basta con pasar el resultado de
# aedv_prediccion_stl_arima(), sin más.
#
# El argumento "tendencia" admite además, explícitamente, el tibble que devuelve
# aedv_tendencia_stl(). Hace falta cuando el resultado ya no trae el atributo
# -un verbo de dplyr lo descarta- o cuando se quiere otra: en un cuadro de mando
# es la componente que de verdad cambia al mover las ventanas del suavizado -la
# suma de tendencia y estacionalidad reconstruye la serie casi igual con
# cualquier ventana-, así que sin ella los controles parecen no hacer nada.
#
# El ggplot de partida se construye SIN la línea naranja (tendencia = NULL):
# aquí se añade como traza de plotly, sobre el eje horizontal ya convertido a
# texto ISO y con su valor incorporado al cuadro emergente. Dibujarla en las dos
# capas la duplicaría.
aedv_dibujar_stl_arima_dinamico <- function(resultado, titulo = NULL,
                                            x_lab = "", y_lab = NULL,
                                            tendencia = .aedv_tendencia_por_omision(resultado)) {
  p <- aedv_dibujar_stl_arima(resultado, x_lab = x_lab, y_lab = y_lab,
                              tendencia = NULL)
  idx <- tsibble::index_var(resultado)
  val <- setdiff(names(resultado), c(idx, "prediccion", "lower", "upper"))[1]
  .aedv_prediccion_a_plotly(p, titulo, .aedv_indice_es_fecha(resultado),
                             resultado, idx, val, "Descomposición",
                             tendencia = tendencia)
}


# ¿El índice temporal de la serie es una fecha de calendario (yearmonth,
# yearquarter, yearweek, Date, POSIXct) o un simple número (un año suelto,
# como en el PIB o los nacimientos)? Lo necesita .aedv_prediccion_a_plotly()
# para saber si merece la pena convertir el eje x a fecha.
.aedv_indice_es_fecha <- function(resultado) {
  idx <- tsibble::index_var(resultado)
  inherits(resultado[[idx]], c("yearmonth", "yearquarter", "yearweek", "Date", "POSIXct"))
}


# Dibuja un mapa coroplético con leaflet a partir de un objeto espacial, un vector
# de valores para colorear por cuantiles, las etiquetas de los tooltips y el título
# de la leyenda. Los argumentos lng, lat y zoom fijan el encuadre inicial del mapa
# (por defecto, una vista del mundo entero).
aedv_mapa_coropletico <- function(geoj,value,region_labels,legend_title,
                                  lng = 25, lat = 22, zoom = 2){
pal <- colorQuantile("YlOrRd", value, n = 9)
p <-  geoj |>
    leaflet() |>
    setView(lng = lng, lat = lat, zoom = zoom)  |>
    addPolygons(
      fillColor = ~pal(value),
      weight = 2,
      opacity = 1,
      color = "white",
      dashArray = "3",
      fillOpacity = 0.7,
      highlightOptions = highlightOptions(
        weight = 2,
        color = rgb(0.2,0.2,0.2),
        dashArray = "",
        fillOpacity = 0.7,
        bringToFront = TRUE
      ),
      label = region_labels
    ) |>
    addLegend("bottomleft",
      pal = pal,
      values = value,
      title = legend_title,
      labFormat = function(type, cuts, p) { # rango de valores de cada tramo
        n = length(cuts)
        x = prettyNum(round(cuts,
            digits=max(5-nchar(as.character(round(max(na.omit(value))))),0)),
            big.mark = ".", decimal.mark = ","
        )
        paste(x[-n], "–", x[-1])
      },
      opacity = 1
    )
  return(p)
}


# ==========================================================================
# Proyecto personal: incorporación de los datasets a la memoria
# ==========================================================================
#
# Cada *_DatasetSelection.Rmd del proyecto personal guarda, junto al rds y al
# xlsx de metadatos del dataset elegido, un fichero <nombre>_analisis.Rmd con
# su sección «Análisis del dataset» preparada para la memoria: el
# identificador del dataset ya fijado y los títulos al nivel de la sección
# «Comprensión de los datos». Ese fichero lo genera aedv_generar_analisis()
# al compilar el *_DatasetSelection.Rmd, y aedv_incluir_dataset() lo inserta
# en la memoria y devuelve los datos.
#
# Marcas que estas dos funciones buscan en los *_DatasetSelection.Rmd (como
# líneas completas, sin nada más en la línea):
#   # Análisis del dataset                                  inicio del análisis
#   <!-- FIN DE LA PARTE QUE SE INCORPORA A LA MEMORIA -->  fin del análisis
#   <!-- INICIO DE LA VALIDACIÓN AEDV -->                   lo que omite
#   <!-- FIN DE LA VALIDACIÓN AEDV -->                      validacion = FALSE


# Genera el fichero <nombre>_analisis.Rmd que aedv_incluir_dataset() inserta
# en la memoria y escribe las instrucciones para el alumno. Se llama desde el
# propio *_DatasetSelection.Rmd, en un chunk con results='asis' situado
# después de la marca de fin del análisis.
#
# Recorta del documento que se está compilando la sección «Análisis del
# dataset» y la adapta para la memoria: aplica las sustituciones (p. ej. fija
# el identificador del dataset), baja los títulos dos niveles para que
# cuelguen de «Comprensión de los datos» (## en la memoria) y pone como título
# de la sección el nombre del dataset.
#
# Argumentos:
#   nombre        : nombre base de los ficheros del dataset (el `filename` del
#                   documento). El fichero generado es <nombre>_analisis.Rmd.
#   titulo        : título de la sección del dataset en la memoria.
#   ficheros      : ficheros de datos que el alumno copia junto a la memoria
#                   (rds y xlsx de metadatos).
#   tablas        : tablas que devuelve el análisis (su `tablas_exportadas`).
#   sugerencia    : nombre que se propone al alumno para el resultado de
#                   aedv_incluir_dataset().
#   desempaquetar : solo con varias tablas. Vector con nombre: nombre que se
#                   propone para cada tabla -> tabla de la lista, p. ej.
#                   c(clima_obs = "data", clima_proy = "proy").
#   comentarios   : comentario de cada línea de desempaquetar (opcional).
#   sustituciones : vector con nombre: patrón (expresión regular) de una línea
#                   del análisis -> línea que la sustituye. Cada patrón tiene
#                   que aparecer exactamente una vez.
aedv_generar_analisis <- function(nombre, titulo, ficheros, tablas = "data",
                                  sugerencia = paste0("datos_", nombre),
                                  desempaquetar = NULL, comentarios = NULL,
                                  sustituciones = NULL) {
  origen  <- knitr::current_input(dir = TRUE)
  carpeta <- basename(dirname(origen))
  fuente  <- sub("\r$", "", readLines(origen, encoding = "UTF-8", warn = FALSE))

  ini <- which(fuente == "# Análisis del dataset")
  fin <- which(fuente == "<!-- FIN DE LA PARTE QUE SE INCORPORA A LA MEMORIA -->")
  if (length(ini) != 1 || length(fin) != 1 || ini > fin)
    stop("No se encuentra en ", basename(origen), " la sección '# Análisis del ",
         "dataset' o su marca de fin", call. = FALSE)
  analisis <- fuente[(ini + 1):(fin - 1)]

  for (patron in names(sustituciones)) {
    pos <- grep(patron, analisis)
    if (length(pos) != 1)
      stop("La línea '", patron, "' aparece ", length(pos), " veces en el ",
           "análisis (tiene que aparecer una)", call. = FALSE)
    analisis[pos] <- sustituciones[[patron]]
  }

  # títulos dos niveles más abajo. Se saltan los chunks y los comentarios
  # HTML, cuyas líneas también pueden empezar por #
  en_chunk <- FALSE
  en_comentario <- FALSE
  for (i in seq_along(analisis)) {
    if (grepl("^\\s*```", analisis[i])) { en_chunk <- !en_chunk; next }
    if (en_chunk) next
    if (grepl("<!--", analisis[i], fixed = TRUE)) en_comentario <- TRUE
    if (!en_comentario && grepl("^#{1,4} ", analisis[i]))
      analisis[i] <- paste0("##", analisis[i])
    if (grepl("-->", analisis[i], fixed = TRUE)) en_comentario <- FALSE
  }

  faltan <- ficheros[!file.exists(ficheros)]
  if (length(faltan) > 0)
    stop("No existen los ficheros: ", paste(faltan, collapse = ", "), call. = FALSE)

  # código que el alumno pega en la memoria; con varias tablas, la lista se
  # separa en una tabla por línea, con su propio nombre
  fichero <- paste0(nombre, "_analisis.Rmd")
  codigo  <- paste0(sugerencia, " <- aedv_incluir_dataset(\"", fichero, "\")")
  if (length(tablas) > 1) {
    if (!setequal(desempaquetar, tablas))
      stop("desempaquetar tiene que nombrar las tablas: ",
           paste(tablas, collapse = ", "), call. = FALSE)
    nombres <- formatC(names(desempaquetar), width = -max(nchar(names(desempaquetar))))
    lineas  <- paste0(nombres, " <- ", sugerencia, "$", desempaquetar)
    if (!is.null(comentarios))
      lineas <- paste0(formatC(lineas, width = -max(nchar(lineas))), "  # ",
                       comentarios[names(desempaquetar)])
    codigo <- c(codigo, paste0("# este dataset tiene ", length(tablas),
                               " tablas: aedv_incluir_dataset() las devuelve en una lista"),
                lineas)
  }

  writeLines(c(
    "<!--",
    paste0("Análisis del dataset \"", titulo, "\", generado por ",
           basename(origen), " el ", Sys.Date(), "."),
    "Se incluye en la sección \"Comprensión de los datos\" de la memoria con este código, en un",
    "chunk con las opciones results='asis' y cache=FALSE:",
    "",
    paste0("    ", codigo),
    "",
    "No conviene editar este fichero: el texto explicativo del dataset se escribe en la memoria,",
    "antes o después del chunk que lo incluye.",
    "-->",
    "",
    paste0("### ", titulo),
    analisis
  ), fichero, useBytes = TRUE)

  # instrucciones para el alumno, con los nombres de sus ficheros
  lista <- function(x) {
    x <- paste0("`", x, "`")
    if (length(x) == 1) x else paste(paste(x[-length(x)], collapse = ", "), "y", x[length(x)])
  }
  resultado <- if (length(tablas) == 1)
    paste0("los datos quedan en la tabla `", sugerencia, "` (un tibble), con la que se ",
           "sigue trabajando en el resto de la memoria.")
  else
    paste0("los datos quedan en las tablas ", lista(names(desempaquetar)), " (tibbles), con ",
           "las que se sigue trabajando en el resto de la memoria. Este dataset tiene ",
           length(tablas), " tablas: `aedv_incluir_dataset()` las devuelve juntas en una ",
           "lista, y las líneas que siguen a la llamada las separan, cada una con su nombre.")
  validacion <- any(analisis == "<!-- INICIO DE LA VALIDACIÓN AEDV -->")

  cat(
    "Al compilar este documento se ha generado el fichero **`", fichero, "`**: es la ",
    "sección *Análisis del dataset* de arriba, preparada para la sección **Comprensión de ",
    "los datos** de la memoria, con los títulos al nivel de esa sección. Para incorporar ",
    "el dataset a la memoria:\n\n",
    "1. Copia a la carpeta de la memoria estos ficheros de la carpeta `", carpeta, "`: ",
    lista(c(ficheros, fichero)), ".\n",
    "2. Comprueba que en la carpeta de la memoria está también `utilidades.R` (viene en la ",
    "carpeta `ProyectoPersonal`, junto al modelo de memoria) y que la memoria lo carga en ",
    "el chunk de librerías: `source(\"utilidades.R\")`.\n",
    "3. Pega este chunk en la sección *Comprensión de los datos* de la memoria, cambiando ",
    "los nombres de las tablas por otros que describan el dataset:\n\n",
    "````\n",
    "```{r, results='asis', cache=FALSE}\n",
    paste0(codigo, "\n"),
    "```\n",
    "````\n\n",
    "Al compilar la memoria, el análisis aparece en el lugar del chunk y ", resultado,
    " Si incorporas varios datasets, pega un chunk por cada uno, con nombres de tabla distintos.",
    if (validacion)
      paste0(" En un dataset auxiliar, que no tiene que cumplir los requisitos de AEDV, ",
             "`aedv_incluir_dataset(\"", fichero, "\", validacion = FALSE)` omite la validación."),
    "\n\n",
    "El texto explicativo del dataset (qué mide, qué significa cada variable y cada una de ",
    "sus opciones...) se escribe en la memoria, antes o después del chunk, no en `", fichero,
    "`: ese fichero se sustituye cada vez que se vuelve a generar.\n",
    sep = ""
  )
  invisible(fichero)
}

# Inserta en la memoria el análisis de un dataset (el fichero
# <nombre>_analisis.Rmd que genera su *_DatasetSelection.Rmd) y devuelve sus
# datos. Se llama desde un chunk con results='asis', sin el cual el análisis
# saldría como texto literal:
#
#   {r, results='asis', cache=FALSE}
#   paro <- aedv_incluir_dataset("35176_analisis.Rmd")
#
# Devuelve un tibble con los datos (la tabla `data` del análisis). Los
# datasets con varias tablas (EDGAR, CLIMA, INFORME PISA) devuelven una lista
# con nombre, una entrada por tabla (p. ej. clima$data y clima$proy). Qué
# tablas devuelve cada análisis lo declara el propio fichero, en la variable
# `tablas_exportadas` de su primer chunk.
#
# El análisis se ejecuta en un entorno propio: los objetos auxiliares que crea
# (tablas intermedias, funciones de la validación...) no pasan a la memoria,
# así que no pueden pisar los del alumno. Y como el resultado se asigna al
# nombre que elija el alumno, incluir varios datasets no sobrescribe ninguno.
#
# Argumentos:
#   fichero    : ruta del fichero <nombre>_analisis.Rmd. Su rds y su xlsx de
#                metadatos tienen que estar en la carpeta de la memoria.
#   validacion : FALSE omite la validación de requisitos AEDV, para los
#                datasets auxiliares que no tienen que cumplirlos.
#
# Fuera de la compilación (al ejecutar el chunk en RStudio con Ctrl+Alt+C) no
# se genera el informe: solo se ejecuta su código para devolver los datos, de
# modo que se puede seguir trabajando chunk a chunk.
aedv_incluir_dataset <- function(fichero, validacion = TRUE) {
  if (!file.exists(fichero))
    stop("No se encuentra el fichero '", fichero, "'. Cópialo a la carpeta de ",
         "la memoria junto con el .rds y el _metadatos.xlsx del mismo dataset.",
         call. = FALSE)

  lineas <- readLines(fichero, encoding = "UTF-8", warn = FALSE)

  # sin validación: se quitan las líneas entre las dos marcas que la delimitan
  if (!validacion) {
    ini <- which(trimws(lineas) == "<!-- INICIO DE LA VALIDACIÓN AEDV -->")
    fin <- which(trimws(lineas) == "<!-- FIN DE LA VALIDACIÓN AEDV -->")
    if (length(ini) != 1 || length(fin) != 1 || ini > fin)
      stop("El fichero '", fichero, "' no tiene marcada la parte de validación: ",
           "no se puede usar validacion = FALSE con él.", call. = FALSE)
    lineas <- lineas[-(ini:fin)]
  }

  # entorno propio, hijo del de la memoria: el análisis ve las librerías y los
  # ajustes de la memoria (p. ej. select <- dplyr::select), pero lo que crea
  # se queda dentro
  entorno <- new.env(parent = parent.frame())
  compilando <- isTRUE(getOption("knitr.in.progress"))

  if (compilando) {
    if (!identical(knitr::opts_current$get("results"), "asis"))
      stop("El chunk que llama a aedv_incluir_dataset() necesita la opción ",
           "results='asis', así: {r, results='asis', cache=FALSE}", call. = FALSE)
    # echo y cache se fijan aquí para que el análisis salga igual aunque la
    # memoria cambie esas opciones por defecto
    informe <- knitr::knit_child(text = lineas, envir = entorno, quiet = TRUE,
                                 options = list(echo = FALSE, cache = FALSE))
    cat(informe, sep = "\n")
  } else {
    # ejecución chunk a chunk: solo el código, con su salida descartada
    codigo <- knitr::purl(text = lineas, quiet = TRUE, documentation = 0)
    utils::capture.output(eval(parse(text = codigo), envir = entorno))
  }

  tablas <- get0("tablas_exportadas", envir = entorno, inherits = FALSE,
                 ifnotfound = "data")
  faltan <- tablas[!vapply(tablas, exists, logical(1), envir = entorno,
                           inherits = FALSE)]
  if (length(faltan) > 0)
    stop("El análisis de '", fichero, "' no ha creado la(s) tabla(s): ",
         paste(faltan, collapse = ", "), call. = FALSE)

  if (length(tablas) == 1)
    return(invisible(get(tablas, envir = entorno)))

  if (!compilando)
    message("Este dataset tiene ", length(tablas), " tablas: se devuelven en ",
            "una lista con ", paste0("$", tablas, collapse = ", "))
  invisible(mget(tablas, envir = entorno))
}
