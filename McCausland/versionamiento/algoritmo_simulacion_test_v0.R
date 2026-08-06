source("algoritmo_simulacion.R")

set.seed(123)

# ============================================================
# Caso 1: m = 4, p = 10  (dimensiones pequeñas)
# ============================================================
sim1 <- modelo_de_simulacion(m = 4, p = 10, n = 1000)

cat("=== Caso 1: m=4, p=10 ===\n")
cat("Dimensión de Y (observaciones):", nrow(sim1$Y), "×", ncol(sim1$Y), "\n")
cat("Dimensión de alpha (estados):  ", nrow(sim1$alpha), "×", ncol(sim1$alpha), "\n")
cat("\nMatriz Q (covarianza de innovaciones):\n")
print(round(sim1$Q, 4))

# ============================================================
# Caso 2: m = 10, p = 100  (data-rich environment)
# ============================================================
sim2 <- modelo_de_simulacion(m = 10, p = 100, n = 1000)

cat("\n=== Caso 2: m=10, p=100 ===\n")
cat("Dimensión de Y:", nrow(sim2$Y), "×", ncol(sim2$Y), "\n")
cat("Dimensión de alpha:", nrow(sim2$alpha), "×", ncol(sim2$alpha), "\n")
cat("\nMatriz Q (covarianza de innovaciones):\n")
print(round(sim2$Q, 4))