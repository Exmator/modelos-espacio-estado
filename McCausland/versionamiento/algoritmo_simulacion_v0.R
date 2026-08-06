# ============================================================
# Algoritmo tomado de:
# McCausland et al. (2011) - Simulation smoothing for state–space models
# Página 205: Computational experiments
# ============================================================

### Algoritmo de generación de estados y observaciones simulados ###
modelo_de_simulacion <- function(m, p, n = 1000) {
  # --- Parametros de entrada ---
  #   m = cantidad de factores del estado
  #   p = cantidad de series observadas
  #   n = número de períodos de tiempo
  
  # --- Parámetros del modelo ---
  iota_m <- matrix(1, nrow = m, ncol = 1)     # vector de unos m×1
  Im     <- diag(m)                           # identidad m×m
  Ip     <- diag(p)                           # identidad p×p
  
  # Matriz de covarianza de innovaciones del estado
  # Q = (0.2)^2 * (1/2 * Im + 1/2 * iota_m * iota_m^T)
  Q  <- (0.2)^2 * (0.5 * Im + 0.5 * iota_m %*% t(iota_m))
  
  Q1 <- Im          # covarianza del estado inicial α₁
  T  <- 0.9 * Im    # matriz de transición
  D  <- Ip          # covarianza del error de medición (diagonal)
  a  <- matrix(0, nrow = m, ncol = 1)  # media inicial
  
  # Factor loading matrix Z (p × m)
  # Zij ~ N(0, (0.001)^2)
  Z <- matrix(rnorm(p * m, mean = 0, sd = 0.001), nrow = p, ncol = m)
  
  # --- Cholesky para simular desde las distribuciones ---
  chol_Q  <- chol(Q)
  chol_Q1 <- chol(Q1)
  chol_D  <- chol(D)
  
  # --- Simular estados αt (m × n) ---
  alpha <- matrix(0, nrow = m, ncol = n)
  
  # Estado inicial: α₁ ~ N(a, Q1)
  alpha[, 1] <- a + t(chol_Q1) %*% rnorm(m)
  
  # Estados siguientes: αt+1 = T * αt + vt,  vt ~ N(0, Q)
  for (t in 1:(n - 1)) {
    v_t         <- t(chol_Q) %*% rnorm(m)
    alpha[, t+1] <- T %*% alpha[, t] + v_t
  }
  
  # --- Simular observaciones yt (p × n) ---
  # yt = Z * αt + ut,  ut ~ N(0, D)
  Y <- matrix(0, nrow = p, ncol = n)
  for (t in 1:n) {
    u_t     <- t(chol_D) %*% rnorm(p)
    Y[, t]  <- Z %*% alpha[, t] + u_t
  }
  
  list(Y = Y, alpha = alpha, Z = Z, T = T, Q = Q, Q1 = Q1, D = D, a = a, m = m, p = p, n = n)
}