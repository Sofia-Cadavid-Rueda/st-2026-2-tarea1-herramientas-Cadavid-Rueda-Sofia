# 03-evaluacion.R
# Funciones de medidas de error y validacion de residuos:

# medidas(e, y)-> MSE, MAD, MAPE, MASE

medidas <- function(e, y, y_entrenamiento = y, s = 1) {
  stopifnot(length(e) == length(y))
  stopifnot(!anyNA(e), !anyNA(y), !anyNA(y_entrenamiento))
  stopifnot(length(y_entrenamiento) > s)
  
  N    <- length(e)
  mse  <- mean(e^2)
  mad  <- mean(abs(e))
  mape <- 100 * mean(abs(e / y))       # y pareado con e: correcto para MAPE
  
  # Denominador del MASE: SIEMPRE sobre y_entrenamiento (el tramo de
  # estimacion), no sobre y (que puede ser el tramo de validacion y tener
  # una longitud distinta a la de e).
  denom_mase <- mean(abs(diff(y_entrenamiento, lag = s)))
  mase <- mad / denom_mase
  
  list(MSE = mse, MAD = mad, MAPE = mape, MASE = mase, N = N)
}

# ljung_box(r, T, m, p)-> Q_m ~ chi2_(m-p)

ljung_box <- function(r, T, m, p = 0) {
  stopifnot(length(r) >= m, m >= 1)
  stopifnot(p >= 0, p < m)
  
  r  <- r[1:m]
  Q  <- T * (T + 2) * sum(r^2 / (T - (1:m)))
  gl <- m - p
  
  valor_critico <- stats::qchisq(0.95, gl)
  valor_p       <- 1 - stats::pchisq(Q, gl)
  
  list(
    estadistico   = Q,
    gl            = gl,
    valor_critico = valor_critico,
    valor_p       = valor_p,
    rechaza_H0    = Q > valor_critico
  )
}

# jarque_bera(e)-> JB ~ chi2_2

jarque_bera <- function(e) {
  stopifnot(!anyNA(e))
  N  <- length(e)
  mu <- mean(e)
  m2 <- mean((e - mu)^2)
  m3 <- mean((e - mu)^3)
  m4 <- mean((e - mu)^4)
  
  A <- m3 / m2^1.5   # asimetria
  K <- m4 / m2^2      # curtosis
  
  JB <- (N / 6) * (A^2 + (K - 3)^2 / 4)
  gl <- 2
  
  valor_critico <- stats::qchisq(0.95, gl)
  valor_p       <- 1 - stats::pchisq(JB, gl)
  
  list(
    estadistico   = JB,
    gl            = gl,
    valor_critico = valor_critico,
    valor_p       = valor_p,
    asimetria     = A,
    curtosis      = K,
    N             = N,
    advertencia_muestra_pequena = N < 20
  )
}

# durbin_watson(e)-> d = sum((e_t - e_{t-1})^2) / sum(e_t^2)

durbin_watson <- function(e) {
  stopifnot(!anyNA(e), length(e) >= 2)
  d <- sum(diff(e)^2) / sum(e^2)
  list(estadistico = d)
}

# validar_errores(e, ...) -> grafico de errores, correlograma, prueba t de
# media cero, Ljung-Box, Jarque-Bera, Durbin-Watson

validar_errores <- function(e, m = NULL, p = 0) {
  stopifnot(!anyNA(e))
  N <- length(e)
  if (is.null(m)) m <- min(floor(N / 4), 24)
  
  grafico_tiempo <- ggplot2::ggplot(
    tibble::tibble(t = seq_along(e), e = e),
    ggplot2::aes(x = t, y = e)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 0, color = "red") +
    ggplot2::labs(x = "t", y = "error", title = "Errores de un paso en el tiempo") +
    ggplot2::theme_minimal()
  
  cg <- correlograma(e, m = m)
  
  prueba_t <- stats::t.test(e, mu = 0)
  lb <- ljung_box(cg$r_h, T = N, m = m, p = p)
  jb <- jarque_bera(e)
  dw <- durbin_watson(e)
  
  list(
    grafico_tiempo      = grafico_tiempo,
    correlograma        = cg,
    prueba_t_media_cero = prueba_t,
    ljung_box           = lb,
    jarque_bera         = jb,
    durbin_watson       = dw,
    N                   = N
  )
}