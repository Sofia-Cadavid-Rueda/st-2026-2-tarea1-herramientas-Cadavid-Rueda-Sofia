# 02-metodos.R
# Los ocho metodos de pronostico, uno por funcion. Cada funcion recibe un
# vector numerico y (ordenado en el tiempo, sin huecos) y devuelve una lista
# con al menos:
#   - yhat: vector misma longitud que y, yhat[t] = pronostico de Y_t hecho
#     con informacion hasta t-1 (posiciones de calentamiento = NA)
#   - pronosticar: funcion de un entero h -> vector de h pronosticos
#     extramuestrales Y_(T+1),...,Y_(T+h)
#   - parametros: lista nombrada con los parametros usados

# ajustar_media(y)

ajustar_media <- function(y) {
  stopifnot(is.numeric(y), !anyNA(y))
  Tn <- length(y)
  stopifnot(Tn >= 2)
  
  # yhat[t] = pronostico recursivo de Y_t con info hasta t-1:
  # Yhat_2 = Y_1 ; Yhat_{t+1} = (1/t) * sum_{i<=t} Y_i
  yhat <- rep(NA_real_, Tn)
  for (t in 2:Tn) {
    yhat[t] <- mean(y[1:(t - 1)])
  }
  
  media_total <- mean(y)  # Ybar sobre toda la muestra, usada para el pronostico extramuestral
  
  pronosticar <- function(h) {
    stopifnot(is.numeric(h), length(h) == 1, h >= 1, h == as.integer(h))
    rep(media_total, h)  # Yhat_{T+h} = Ybar para todo h
  }
  
  list(
    yhat        = yhat,
    pronosticar = pronosticar,
    parametros  = list(media_muestral = media_total)
  )
}

#   ajustar_mm(y, k)

ajustar_mm <- function(y, k) {
  stopifnot(is.numeric(y), !anyNA(y))
  Tn <- length(y)
  stopifnot(is.numeric(k), length(k) == 1, k == as.integer(k))
  stopifnot(k >= 2, k < Tn)
  
  # yhat[t] = MM_{t-1}(k), definida para t-1 >= k, es decir t >= k+1
  yhat <- rep(NA_real_, Tn)
  for (t in (k + 1):Tn) {
    yhat[t] <- mean(y[(t - k):(t - 1)])
  }
  
  mm_final <- mean(y[(Tn - k + 1):Tn])  # MM_T(k), nivel al final de la muestra
  
  pronosticar <- function(h) {
    stopifnot(is.numeric(h), length(h) == 1, h >= 1, h == as.integer(h))
    rep(mm_final, h)  # todo horizonte extramuestral pronostica el mismo numero
  }
  
  list(
    yhat        = yhat,
    pronosticar = pronosticar,
    parametros  = list(k = k, mm_T = mm_final)
  )
}

# ajustar_ses(y, alpha)

ajustar_ses <- function(y, alpha) {
  stopifnot(is.numeric(y), !anyNA(y))
  Tn <- length(y)
  stopifnot(Tn >= 2)
  stopifnot(is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha < 1)
  
  # yhat[t] = pronostico de Y_t con info hasta t-1
  # Yhat_2 = Y_1 ; Yhat_{t+1} = alpha*Y_t + (1-alpha)*Yhat_t
  yhat <- rep(NA_real_, Tn)
  yhat[2] <- y[1]
  for (t in 2:(Tn - 1)) {
    yhat[t + 1] <- alpha * y[t] + (1 - alpha) * yhat[t]
  }
  
  # Verificacion: la forma de correccion de error coincide con la de
  # promedio ponderado (pesos geometricos + peso residual sobre Y_1,
  # que viene de la inicializacion Yhat_2 = Y_1)
  t_verif <- Tn - 1
  i <- 0:(t_verif - 2)
  suma_pesos <- alpha * sum((1 - alpha)^i * y[t_verif - i])
  yhat_pesos <- suma_pesos + (1 - alpha)^(t_verif - 1) * y[1]
  dif_max <- abs(yhat_pesos - yhat[t_verif + 1])
  stopifnot(dif_max < 1e-9)
  
  pronostico_final <- alpha * y[Tn] + (1 - alpha) * yhat[Tn]
  
  pronosticar <- function(h) {
    stopifnot(is.numeric(h), length(h) == 1, h >= 1, h == as.integer(h))
    rep(pronostico_final, h)  # constante en h, igual que media y MM
  }
  
  list(
    yhat        = yhat,
    pronosticar = pronosticar,
    parametros  = list(alpha = alpha, yhat_T1 = pronostico_final,
                       dif_max_verificacion = dif_max)
  )
}

# ajustar_dmm(y, k)

ajustar_dmm <- function(y, k) {
  stopifnot(is.numeric(y), !anyNA(y))
  Tn <- length(y)
  stopifnot(is.numeric(k), length(k) == 1, k == as.integer(k))
  stopifnot(k >= 2, 2 * k <= Tn)
  
  # MM_t(k), definida para t = k,...,Tn
  mm <- rep(NA_real_, Tn)
  for (t in k:Tn) {
    mm[t] <- mean(y[(t - k + 1):t])
  }
  
  # DMM_t(k): media movil de la propia MM, definida para t = 2k-1,...,Tn
  dmm <- rep(NA_real_, Tn)
  for (t in (2 * k - 1):Tn) {
    dmm[t] <- mean(mm[(t - k + 1):t])
  }
  
  # E_t (nivel corregido) y beta1_t (pendiente local)
  idx <- (2 * k - 1):Tn
  E_t <- rep(NA_real_, Tn)
  beta1_t <- rep(NA_real_, Tn)
  E_t[idx]     <- 2 * mm[idx] - dmm[idx]
  beta1_t[idx] <- (2 / (k - 1)) * (mm[idx] - dmm[idx])
  
  # yhat[t] = pronostico de Y_t con info hasta t-1: Yhat_{t+1} = E_t + beta1_t
  # definido desde t = 2k
  yhat <- rep(NA_real_, Tn)
  for (t in (2 * k):Tn) {
    yhat[t] <- E_t[t - 1] + beta1_t[t - 1]
  }
  
  E_T     <- E_t[Tn]
  beta1_T <- beta1_t[Tn]
  
  pronosticar <- function(h) {
    stopifnot(is.numeric(h), length(h) == 1, h >= 1, h == as.integer(h))
    E_T + beta1_T * (1:h)  # extrapolacion lineal, depende de h
  }
  
  list(
    yhat        = yhat,
    pronosticar = pronosticar,
    parametros  = list(k = k, E_T = E_T, beta1_T = beta1_T,
                       MM_t = mm, DMM_t = dmm, E_t = E_t, beta1_t = beta1_t)
  )
}

#   ajustar_tendencia(y, tipo)      # tipo: "lineal" | "cuadratica" | "exponencial"

ajustar_tendencia <- function(y, tipo = c("lineal", "cuadratica", "exponencial"),
                              corregir_sesgo = FALSE) {
  tipo <- match.arg(tipo)
  stopifnot(is.numeric(y), !anyNA(y))
  Tn <- length(y)
  t_idx <- seq_len(Tn)
  
  if (tipo == "exponencial" && any(y <= 0)) {
    stop("ajustar_tendencia: tipo 'exponencial' requiere y > 0 en todas las posiciones")
  }
  
  # variable respuesta de la regresion: log(y) si es exponencial, y si no
  y_reg <- if (tipo == "exponencial") log(y) else y
  
  X <- switch(tipo,
              lineal      = cbind(Intercepto = 1, t = t_idx),
              cuadratica  = cbind(Intercepto = 1, t = t_idx, t2 = t_idx^2),
              exponencial = cbind(Intercepto = 1, t = t_idx)   # a + theta*t
  )
  k <- ncol(X)
  stopifnot(Tn > k)
  
  # --- Ecuaciones normales
  XtX  <- crossprod(X)
  Xty  <- crossprod(X, y_reg)
  beta <- as.numeric(solve(XtX, Xty))
  names(beta) <- colnames(X)
  
  ajustado_reg <- as.numeric(X %*% beta)   # ajustado en escala de la regresion
  e_reg <- y_reg - ajustado_reg            # residuos en escala de la regresion
  gl   <- Tn - k
  sigma2 <- sum(e_reg^2) / gl
  
  # yhat en la escala ORIGINAL de y, para que sea comparable con los demas metodos
  factor_sesgo <- if (tipo == "exponencial" && corregir_sesgo) exp(sigma2 / 2) else 1
  yhat <- if (tipo == "exponencial") exp(ajustado_reg) * factor_sesgo else ajustado_reg
  e    <- y - yhat   # error en escala original, para medidas() y validar_errores()
  
  XtX_inv <- solve(XtX)
  se_ols  <- sqrt(diag(sigma2 * XtX_inv))
  
  # --- Error estandar robusto HAC (Newey-West, nucleo de Bartlett), sobre e_reg ---
  L  <- floor(4 * (Tn / 100)^(2 / 9))
  Xe <- X * e_reg
  Meat <- crossprod(Xe)
  if (L >= 1) {
    for (l in 1:L) {
      w <- 1 - l / (L + 1)
      Gamma_l <- crossprod(Xe[(l + 1):Tn, , drop = FALSE],
                           Xe[1:(Tn - l), , drop = FALSE])
      Meat <- Meat + w * (Gamma_l + t(Gamma_l))
    }
  }
  Var_hac <- XtX_inv %*% Meat %*% XtX_inv
  se_hac  <- sqrt(diag(Var_hac))
  t_hac <- beta / se_hac
  p_hac <- 2 * (1 - stats::pt(abs(t_hac), gl))
  
  tabla_coef <- tibble::tibble(
    termino     = names(beta),
    estimacion  = beta,
    ee_ols      = se_ols,
    ee_hac      = se_hac,
    t_hac       = t_hac,
    valor_p_hac = p_hac
  )
  
  SCT <- sum((y_reg - mean(y_reg))^2)
  R2  <- 1 - sum(e_reg^2) / SCT
  dw  <- sum(diff(e_reg)^2) / sum(e_reg^2)
  
  # --- Verificacion contra lm() ---
  coef_lm <- unname(stats::coef(stats::lm(y_reg ~ X[, -1, drop = FALSE])))
  dif_max <- max(abs(beta - coef_lm))
  stopifnot(dif_max < 1e-8)
  
  pronosticar <- function(h) {
    stopifnot(is.numeric(h), length(h) == 1, h >= 1, h == as.integer(h))
    tf <- Tn + seq_len(h)
    val_reg <- switch(tipo,
                      lineal      = beta[1] + beta[2] * tf,
                      cuadratica  = beta[1] + beta[2] * tf + beta[3] * tf^2,
                      exponencial = beta[1] + beta[2] * tf
    )
    if (tipo == "exponencial") exp(val_reg) * factor_sesgo else val_reg
  }
  
  list(
    yhat        = yhat,
    pronosticar = pronosticar,
    parametros  = list(
      tipo = tipo, coeficientes = beta, tabla = tabla_coef,
      beta0 = if (tipo == "exponencial") exp(beta[1]) else NA,
      beta1 = if (tipo == "exponencial") exp(beta[2]) else NA,
      corregir_sesgo = corregir_sesgo, factor_sesgo = factor_sesgo,
      R2 = R2, sigma2 = sigma2, dw = dw, L_hac = L, gl = gl,
      dif_max_verificacion = dif_max
    )
  )
}

# ajustar_holt(y, alpha, beta)

ajustar_holt <- function(y, alpha, beta) {
  stopifnot(is.numeric(y), !anyNA(y))
  Tn <- length(y)
  stopifnot(Tn >= 3)
  stopifnot(is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha < 1)
  stopifnot(is.numeric(beta),  length(beta)  == 1, beta  > 0, beta  < 1)
  
  L  <- rep(NA_real_, Tn)
  Th <- rep(NA_real_, Tn)   # T-hat
  
  # Calentamiento (Def. de las notas): L_1 = Y_1, That_1 = 0
  L[1]  <- y[1]
  Th[1] <- 0
  
  # yhat[t] = pronostico de Y_t con info hasta t-1: Yhat_{t+1} = L_t + That_t
  yhat <- rep(NA_real_, Tn)
  
  for (t in 2:Tn) {
    L[t]  <- alpha * y[t] + (1 - alpha) * (L[t - 1] + Th[t - 1])
    Th[t] <- beta * (L[t] - L[t - 1]) + (1 - beta) * Th[t - 1]
    if (t < Tn) yhat[t + 1] <- L[t] + Th[t]
  }
  
  # Verificacion: forma de correccion de error
  # L_t = L_{t-1} + That_{t-1} + alpha*e_t ; That_t = That_{t-1} + alpha*beta*e_t
  e_verif <- y[2:Tn] - yhat[2:Tn]     # e_t para t = 2,...,Tn (NA en la ultima pos si aplica)
  ok <- which(!is.na(e_verif))
  t_idx <- (2:Tn)[ok]
  e_t   <- e_verif[ok]
  
  L_check  <- L[t_idx - 1] + Th[t_idx - 1] + alpha * e_t
  Th_check <- Th[t_idx - 1] + alpha * beta * e_t
  
  dif_max_L  <- max(abs(L_check  - L[t_idx]))
  dif_max_Th <- max(abs(Th_check - Th[t_idx]))
  stopifnot(dif_max_L < 1e-9, dif_max_Th < 1e-9)
  
  L_T  <- L[Tn]
  Th_T <- Th[Tn]
  
  pronosticar <- function(h) {
    stopifnot(is.numeric(h), length(h) == 1, h >= 1, h == as.integer(h))
    L_T + Th_T * (1:h)   # depende de h, como DMM
  }
  
  list(
    yhat        = yhat,
    pronosticar = pronosticar,
    parametros  = list(
      alpha = alpha, beta = beta, L_T = L_T, Th_T = Th_T,
      L_t = L, That_t = Th,
      dif_max_verificacion = max(dif_max_L, dif_max_Th)
    )
  )
}

# optimizar(y, metodo, rejilla)

optimizar <- function(y, metodo, rejilla) {
  metodo <- match.arg(metodo, c("mm", "dmm", "ses", "holt"))
  stopifnot(is.numeric(y), !anyNA(y))
  
  mse_un_paso <- function(yhat) {
    e  <- y - yhat
    ok <- !is.na(e)
    stopifnot(any(ok))
    mean(e[ok]^2)
  }
  
  if (metodo %in% c("mm", "dmm")) {
    stopifnot(is.numeric(rejilla), length(rejilla) >= 1)
    ajustar_fn <- if (metodo == "mm") ajustar_mm else ajustar_dmm
    
    mse_vec <- vapply(rejilla, function(k) {
      mse_un_paso(ajustar_fn(y, k)$yhat)
    }, numeric(1))
    
    tabla  <- tibble::tibble(k = rejilla, MSE = mse_vec)
    optimo <- tabla[which.min(tabla$MSE), ]
    en_borde <- isTRUE(optimo$k == min(rejilla)) || isTRUE(optimo$k == max(rejilla))
    
  } else if (metodo == "ses") {
    stopifnot(is.numeric(rejilla), length(rejilla) >= 1)
    
    mse_vec <- vapply(rejilla, function(a) {
      mse_un_paso(ajustar_ses(y, a)$yhat)
    }, numeric(1))
    
    tabla  <- tibble::tibble(alpha = rejilla, MSE = mse_vec)
    optimo <- tabla[which.min(tabla$MSE), ]
    en_borde <- isTRUE(optimo$alpha == min(rejilla)) || isTRUE(optimo$alpha == max(rejilla))
    
  } else if (metodo == "holt") {
    stopifnot(is.data.frame(rejilla), all(c("alpha", "beta") %in% names(rejilla)))
    
    mse_vec <- vapply(seq_len(nrow(rejilla)), function(i) {
      mse_un_paso(ajustar_holt(y, rejilla$alpha[i], rejilla$beta[i])$yhat)
    }, numeric(1))
    
    tabla  <- tibble::tibble(alpha = rejilla$alpha, beta = rejilla$beta, MSE = mse_vec)
    optimo <- tabla[which.min(tabla$MSE), ]
    
    rejilla_alpha <- unique(rejilla$alpha)
    rejilla_beta  <- unique(rejilla$beta)
    en_borde_alpha <- isTRUE(optimo$alpha == min(rejilla_alpha)) || isTRUE(optimo$alpha == max(rejilla_alpha))
    en_borde_beta  <- isTRUE(optimo$beta  == min(rejilla_beta))  || isTRUE(optimo$beta  == max(rejilla_beta))
    en_borde <- c(alpha = en_borde_alpha, beta = en_borde_beta)
  }
  
  list(tabla = tabla, optimo = optimo, en_borde = en_borde)
}
