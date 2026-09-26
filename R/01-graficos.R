# 01-graficos.R
# Funcion: graficar_serie(datos, titulo)

graficar_serie <- function(datos, titulo) {
  stopifnot(all(c("fecha", "y") %in% names(datos)))
  stopifnot(is.character(titulo), length(titulo) == 1)
  
  fuente <- attr(datos, "fuente")
  unidad <- attr(datos, "unidad")
  n_obs  <- nrow(datos)
  
  ggplot2::ggplot(datos, ggplot2::aes(x = fecha, y = y)) +
    ggplot2::geom_line(color = "steelblue") +
    ggplot2::scale_x_date(date_labels = "%Y", date_breaks = "2 years") +
    ggplot2::labs(
      title   = titulo,
      x       = "Fecha",
      y       = if (!is.null(unidad)) unidad else "y",
      caption = sprintf(
        "Fuente: %s. n = %d observaciones.",
        if (!is.null(fuente)) fuente else "no declarada",
        n_obs
      )
    ) +
    ggplot2::theme_minimal()
}

# Funcion: correlograma(datos, m)

correlograma <- function(datos_o_y, m = NULL) {
  
  if (is.data.frame(datos_o_y)) {
    stopifnot("y" %in% names(datos_o_y))
    y <- datos_o_y$y
  } else {
    y <- as.numeric(datos_o_y)
  }
  stopifnot(!anyNA(y))
  
  Tn <- length(y)
  if (is.null(m)) m <- min(floor(Tn / 4), 24)
  stopifnot(m >= 1, m < Tn)
  
  ybar  <- mean(y)
  denom <- sum((y - ybar)^2)
  
  r_h <- vapply(1:m, function(h) {
    num <- sum((y[(h + 1):Tn] - ybar) * (y[1:(Tn - h)] - ybar))
    num / denom
  }, numeric(1))
  
  pacf_obj <- stats::pacf(y, lag.max = m, plot = FALSE)
  p_h <- as.numeric(pacf_obj$acf)
  
  banda <- stats::qnorm(0.975) / sqrt(Tn)
  
  dat_acf  <- tibble::tibble(h = 1:m, valor = r_h, tipo = "ACF")
  dat_pacf <- tibble::tibble(h = 1:m, valor = p_h, tipo = "PACF")
  
  plot_uno <- function(dat, etiqueta) {
    ggplot2::ggplot(dat, ggplot2::aes(x = h, y = valor)) +
      ggplot2::geom_hline(yintercept = 0, color = "grey40") +
      ggplot2::geom_hline(yintercept = c(-banda, banda),
                          linetype = "dashed", color = "steelblue") +
      ggplot2::geom_segment(ggplot2::aes(xend = h, yend = 0)) +
      ggplot2::labs(x = "Rezago (h)", y = etiqueta) +
      ggplot2::theme_minimal()
  }
  
  panel <- patchwork::wrap_plots(
    plot_uno(dat_acf, "ACF"),
    plot_uno(dat_pacf, "PACF"),
    ncol = 1
  )
  
  r_h_stats <- as.numeric(stats::acf(y, lag.max = m, plot = FALSE)$acf)[-1]
  dif_max   <- max(abs(r_h - r_h_stats))
  stopifnot(dif_max < 1e-12)
  
  list(
    grafico = panel,
    r_h     = r_h,
    p_h     = p_h,
    banda   = banda,
    m       = m,
    n       = Tn,
    dif_max_verificacion = dif_max
  )
}



