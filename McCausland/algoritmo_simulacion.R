# ============================================================
# Algoritmo tomado de:
# McCausland et al. (2011) - Simulation smoothing for state–space models
# Página 205: Computational experiments
# ============================================================

### Función de validación de dimensiones y tipos de entrada ###
validar_entradas <- function(Q, Q1, T_mat, D, a, Z, m, p, n) {

  # --- Escalares: m, p, n ---
  if (!is.numeric(m) || length(m) != 1 || m != floor(m) || m < 1)
    stop("'m' debe ser un entero positivo.")
  if (!is.numeric(p) || length(p) != 1 || p != floor(p) || p < 1)
    stop("'p' debe ser un entero positivo.")
  if (!is.numeric(n) || length(n) != 1 || n != floor(n) || n < 2)
    stop("'n' debe ser un entero >= 2.")

  # --- Z: matriz p × m ---
  if (!is.matrix(Z))
    stop("'Z' debe ser una matriz.")
  if (!all(dim(Z) == c(p, m)))
    stop(sprintf("'Z' debe ser de dimensión p×m = %d×%d, pero es %d×%d.",
                 p, m, nrow(Z), ncol(Z)))

  # --- T_mat: matriz m × m ---
  if (!is.matrix(T_mat))
    stop("'T_mat' debe ser una matriz.")
  if (!all(dim(T_mat) == c(m, m)))
    stop(sprintf("'T_mat' debe ser de dimensión m×m = %d×%d, pero es %d×%d.",
                 m, m, nrow(T_mat), ncol(T_mat)))

  # --- Q: matriz m × m, simétrica y definida positiva ---
  if (!is.matrix(Q))
    stop("'Q' debe ser una matriz.")
  if (!all(dim(Q) == c(m, m)))
    stop(sprintf("'Q' debe ser de dimensión m×m = %d×%d, pero es %d×%d.",
                 m, m, nrow(Q), ncol(Q)))
  if (!isSymmetric(Q, tol = 1e-10))
    stop("'Q' debe ser una matriz simétrica.")
  if (any(eigen(Q, only.values = TRUE)$values <= 0))
    stop("'Q' debe ser definida positiva (todos sus valores propios > 0).")

  # --- Q1: matriz m × m, simétrica y definida positiva ---
  if (!is.matrix(Q1))
    stop("'Q1' debe ser una matriz.")
  if (!all(dim(Q1) == c(m, m)))
    stop(sprintf("'Q1' debe ser de dimensión m×m = %d×%d, pero es %d×%d.",
                 m, m, nrow(Q1), ncol(Q1)))
  if (!isSymmetric(Q1, tol = 1e-10))
    stop("'Q1' debe ser una matriz simétrica.")
  if (any(eigen(Q1, only.values = TRUE)$values <= 0))
    stop("'Q1' debe ser definida positiva (todos sus valores propios > 0).")

  # --- D: matriz p × p, simétrica y definida positiva ---
  if (!is.matrix(D))
    stop("'D' debe ser una matriz.")
  if (!all(dim(D) == c(p, p)))
    stop(sprintf("'D' debe ser de dimensión p×p = %d×%d, pero es %d×%d.",
                 p, p, nrow(D), ncol(D)))
  if (!isSymmetric(D, tol = 1e-10))
    stop("'D' debe ser una matriz simétrica.")
  if (any(eigen(D, only.values = TRUE)$values <= 0))
    stop("'D' debe ser definida positiva (todos sus valores propios > 0).")

  # --- a: vector m × 1 ---
  if (!is.matrix(a))
    stop("'a' debe ser una matriz columna (m×1).")
  if (!all(dim(a) == c(m, 1)))
    stop(sprintf("'a' debe ser de dimensión m×1 = %d×1, pero es %d×%d.",
                 m, nrow(a), ncol(a)))

  invisible(TRUE)
}

### Algoritmo de generación de estados y observaciones simulados ###
modelo_de_simulacion <- function(Q, Q1, T_mat, D, a, Z, m, p, n = 1000) {
  # --- Parámetros de entrada ---
  #   Q     = covarianza de las innovaciones del estado (m×m)
  #   Q1    = covarianza del estado inicial α₁         (m×m)
  #   T_mat = matriz de transición                     (m×m)
  #   D     = covarianza del error de medición         (p×p)
  #   a     = media del estado inicial                 (m×1)
  #   Z     = factor loading matrix                    (p×m)
  #   m     = cantidad de factores del estado
  #   p     = cantidad de series observadas
  #   n     = número de períodos de tiempo

  # --- Validación de entradas ---
  validar_entradas(Q, Q1, T_mat, D, a, Z, m, p, n)

  # --- Cholesky para simular desde las distribuciones ---
  chol_Q  <- chol(Q)
  chol_Q1 <- chol(Q1)
  chol_D  <- chol(D)

  # --- Simular estados αt (m × n) ---
  alpha       <- matrix(0, nrow = m, ncol = n)
  
  # Estado inicial: α₁ ~ N(a, Q1)
  alpha[, 1]  <- a + t(chol_Q1) %*% rnorm(m)

  # Estados siguientes: αt+1 = T_mat * αt + vt,  vt ~ N(0, Q)
  for (t in 1:(n - 1)) {
    v_t          <- t(chol_Q) %*% rnorm(m)
    alpha[, t+1] <- T_mat %*% alpha[, t] + v_t
  }

  # --- Simular observaciones yt (p × n) ---
  # yt = Z * αt + ut,  ut ~ N(0, D)
  Y <- matrix(0, nrow = p, ncol = n)
  for (t in 1:n) {
    u_t    <- t(chol_D) %*% rnorm(p)
    Y[, t] <- Z %*% alpha[, t] + u_t
  }

  list(Y = Y, alpha = alpha, Z = Z, T_mat = T_mat,
       Q = Q, Q1 = Q1, D = D, a = a, m = m, p = p, n = n)
}
