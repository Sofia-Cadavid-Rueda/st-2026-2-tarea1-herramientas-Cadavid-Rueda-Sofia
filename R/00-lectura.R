# 00-lectura.R
# Funcion: leer_serie(x, fuente, unidad)

leer_serie <- function(x, fuente, unidad) {
  stopifnot(is.character(fuente), length(fuente) == 1)
  stopifnot(is.character(unidad), length(unidad) == 1)
  
  if (inherits(x, "ts")) {
    
    y <- as.numeric(x)
    frecuencia <- stats::frequency(x)
    inicio <- stats::start(x)
    anio_ini <- inicio[1]
    sub_ini  <- if (length(inicio) > 1) inicio[2] else 1
    
    paso_meses <- switch(
      as.character(frecuencia),
      "1"  = 12,   # anual
      "4"  = 3,    # trimestral
      "12" = 1,    # mensual
      stop("leer_serie: frecuencia no soportada: ", frecuencia)
    )
    
    mes_ini <- if (frecuencia == 1) 1 else (sub_ini - 1) * paso_meses + 1
    fecha_ini <- as.Date(sprintf("%d-%02d-01", anio_ini, mes_ini))
    fecha <- seq(fecha_ini, by = paste(paso_meses, "months"), length.out = length(y))
    
  } else if (is.character(x) && length(x) == 1 && file.exists(x)) {
    
    datos_csv <- utils::read.csv(x, stringsAsFactors = FALSE)
    stopifnot(all(c("fecha", "valor") %in% names(datos_csv)))
    
    fecha <- as.Date(datos_csv$fecha)
    y <- as.numeric(datos_csv$valor)
    stopifnot(!anyNA(fecha), !anyNA(y))
    
    d <- as.numeric(diff(fecha))
    paso_tipico <- as.numeric(names(sort(table(d), decreasing = TRUE))[1])
    
    frecuencia <- dplyr::case_when(
      paso_tipico == 1                       ~ 365,
      paso_tipico >= 6   & paso_tipico <= 8   ~ 52,
      paso_tipico >= 27  & paso_tipico <= 31  ~ 12,
      paso_tipico >= 89  & paso_tipico <= 92  ~ 4,
      paso_tipico >= 360 & paso_tipico <= 366 ~ 1,
      TRUE ~ NA_real_
    )
    if (is.na(frecuencia)) {
      stop("leer_serie: no se pudo inferir una frecuencia calendario a partir de las fechas")
    }
    
  } else {
    stop("leer_serie: x debe ser un objeto ts o la ruta a un archivo csv existente ",
         "con columnas 'fecha' y 'valor'")
  }
  
  if (!all(diff(fecha) > 0)) {
    stop("leer_serie: las fechas no son estrictamente crecientes")
  }
  
  if (!inherits(x, "ts")) {
    dias_por_paso <- 365 / frecuencia
    d <- as.numeric(diff(fecha))
    if (any(abs(d - dias_por_paso) > max(3, 0.15 * dias_por_paso))) {
      stop("leer_serie: las fechas no estan equiespaciadas segun la frecuencia inferida")
    }
  }
  
  datos <- tibble::tibble(
    t     = seq_along(y),
    fecha = fecha,
    y     = y
  )
  
  attr(datos, "frecuencia") <- frecuencia
  attr(datos, "fuente")     <- fuente
  attr(datos, "unidad")     <- unidad
  
  datos
}