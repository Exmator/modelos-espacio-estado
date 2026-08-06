source("algoritmo_simulacion.R")

set.seed(123)

# ============================================================
# Parámetros del modelo según McCausland et al. (2011), pág. 205
# Modelo: Dynamic Factor Model
#
# Ecuación de observación: yt = Z*αt + ut,  ut ~ N(0, D)
# Ecuación de estado:   αt+1 = T*αt + vt,  vt ~ N(0, Q)
# Estado inicial:          α1 ~ N(a, Q1)
# ============================================================

# --- Caso 1: m = 4, p = 10 ---
m1 <- 4
p1 <- 10
n  <- 1000

iota_m1 <- matrix(1, nrow = m1, ncol = 1)

# Q = (0.2)^2 * (1/2 * Im + 1/2 * iota_m * iota_m^T)
Q1_caso1  <- (0.2)^2 * (0.5 * diag(m1) + 0.5 * iota_m1 %*% t(iota_m1))

# Q1 = Im  (covarianza del estado inicial)
Q1_ini1   <- diag(m1)

# T = 0.9 * Im  (matriz de transición)
T_mat1    <- 0.9 * diag(m1)

# D = Ip  (covarianza del ruido de medición)
D1        <- diag(p1)

# a = 0_m  (media del estado inicial)
a1        <- matrix(0, nrow = m1, ncol = 1)

# Z: factor loadings, Zij ~ N(0, (0.001)^2)
Z1        <- matrix(rnorm(p1 * m1, mean = 0, sd = 0.001), nrow = p1, ncol = m1)

# ============================================================
# Caso 1: m = 4, p = 10  (dimensiones pequeñas)
# ============================================================
sim1 <- modelo_de_simulacion(
  Q     = Q1_caso1,
  Q1    = Q1_ini1,
  T_mat = T_mat1,
  D     = D1,
  a     = a1,
  Z     = Z1,
  m     = m1,
  p     = p1,
  n     = n
)

cat("=== Caso 1: m=4, p=10 ===\n")
cat("Dimensión de Y (observaciones):", nrow(sim1$Y), "×", ncol(sim1$Y), "\n")
cat("Dimensión de alpha (estados):  ", nrow(sim1$alpha), "×", ncol(sim1$alpha), "\n")
cat("\nMatriz Q (covarianza de innovaciones):\n")
print(round(sim1$Q, 4))

# --- Caso 2: m = 10, p = 100 ---
m2 <- 10
p2 <- 100

iota_m2 <- matrix(1, nrow = m2, ncol = 1)

Q1_caso2  <- (0.2)^2 * (0.5 * diag(m2) + 0.5 * iota_m2 %*% t(iota_m2))
Q1_ini2   <- diag(m2)
T_mat2    <- 0.9 * diag(m2)
D2        <- diag(p2)
a2        <- matrix(0, nrow = m2, ncol = 1)
Z2        <- matrix(rnorm(p2 * m2, mean = 0, sd = 0.001), nrow = p2, ncol = m2)

# ============================================================
# Caso 2: m = 10, p = 100  (data-rich environment)
# ============================================================
sim2 <- modelo_de_simulacion(
  Q     = Q1_caso2,
  Q1    = Q1_ini2,
  T_mat = T_mat2,
  D     = D2,
  a     = a2,
  Z     = Z2,
  m     = m2,
  p     = p2,
  n     = n
)

cat("\n=== Caso 2: m=10, p=100 ===\n")
cat("Dimensión de Y:", nrow(sim2$Y), "×", ncol(sim2$Y), "\n")
cat("Dimensión de alpha:", nrow(sim2$alpha), "×", ncol(sim2$alpha), "\n")
cat("\nMatriz Q (covarianza de innovaciones):\n")
print(round(sim2$Q, 4))

# ============================================================
# Pruebas de validación — entradas incorrectas
# ============================================================
cat("\n=== Pruebas de validación ===\n")

probar_error <- function(descripcion, expr) {
  resultado <- tryCatch(expr, error = function(e) e$message)
  cat(sprintf("  [%s]\n  Error esperado: %s\n\n", descripcion, resultado))
}

# Z con dimensiones incorrectas
probar_error("Z con dimensiones incorrectas",
  modelo_de_simulacion(Q=Q1_caso1, Q1=Q1_ini1, T_mat=T_mat1, D=D1,
                       a=a1, Z=matrix(0, 5, 3), m=m1, p=p1, n=n))

# Q no definida positiva
Q_mala <- Q1_caso1; Q_mala[1,1] <- -1
probar_error("Q no definida positiva",
  modelo_de_simulacion(Q=Q_mala, Q1=Q1_ini1, T_mat=T_mat1, D=D1,
                       a=a1, Z=Z1, m=m1, p=p1, n=n))

# a con dimensión incorrecta
probar_error("a con dimensión incorrecta",
  modelo_de_simulacion(Q=Q1_caso1, Q1=Q1_ini1, T_mat=T_mat1, D=D1,
                       a=matrix(0, m1+1, 1), Z=Z1, m=m1, p=p1, n=n))

# n < 2
probar_error("n < 2",
  modelo_de_simulacion(Q=Q1_caso1, Q1=Q1_ini1, T_mat=T_mat1, D=D1,
                       a=a1, Z=Z1, m=m1, p=p1, n=1))
