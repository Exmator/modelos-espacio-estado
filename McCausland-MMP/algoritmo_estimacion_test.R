# ============================================================
# Test del algoritmo de estimación MMP
# McCausland et al. (2011) - Simulation smoothing for state-space models
#
# Compara estados reales vs estimados con IC al 95% (alpha=0.05),
# usando 1000 draws del algoritmo MMP por individuo. Sobre 100 individuos
# independientes (misma especificación del modelo, distinta
# trayectoria simulada cada vez).
#
# Caso simple: m=2, p=2 para reducir tiempo de cómputo.
# ============================================================

source("McCausland-MMP/algoritmo_simulacion.R")
source("McCausland-MMP/algoritmo_estimacion.R")

set.seed(120)

# ============================================================
# 1. Definir parámetros del modelo (página 205 del artículo)
# ============================================================
m <- 2
p <- 2
n <- 200
 
iota_m <- matrix(1, nrow = m, ncol = 1)

# Q = (0.2)^2 * (1/2 * Im + 1/2 * iota_m * iota_m^T)
Q_sim  <- (0.2)^2 * (0.5 * diag(m) + 0.5 * iota_m %*% t(iota_m))

# Q1 = Im  (covarianza del estado inicial)
Q1_sim <- diag(m)

# T = 0.9 * Im  (matriz de transición)
T_sim  <- 0.9 * diag(m)

# D = Ip  (covarianza del ruido de medición)
D_sim  <- diag(p)

# a = 0_m  (media del estado inicial)
a_sim  <- matrix(0, nrow = m, ncol = 1)

# Z con señal fuerte para que el suavizador pueda recuperar los estados
# (sd=1 en lugar de sd=0.001 del artículo, que era solo para benchmarking)
Z_sim  <- matrix(rnorm(p * m, mean = 0, sd = 1), nrow = p, ncol = m)

cat("=== Parámetros del modelo ===\n")
cat("m (factores):", m, "\n")
cat("p (observaciones):", p, "\n")
cat("n (períodos):", n, "\n")
cat("\nMatriz Z (factor loadings):\n")
print(round(Z_sim, 4))
cat("\nMatriz Q (covarianza innovaciones):\n")
print(round(Q_sim, 4))

# ============================================================
# 2. Función auxiliar: IC empírico al 95% a partir de los draws
#
# Adaptado de Graphcyd.R (ICstates) y ZhuMcStates.R:
# Para cada factor j y cada tiempo t, los Ns draws forman
# una distribución empírica. Se extraen los cuantiles 2.5% y 97.5%
# como límites del IC, y la media empírica como estimación central.
#
# draws[j, t, ] es el vector de Ns valores del factor j en el tiempo t.
# ============================================================
calcular_ic <- function(draws, j, alpha_sig) {
  # Para el factor j, calcula media y límites del IC para cada t
  # draws: array m×n×Ns
  n_t <- dim(draws)[2]

  media <- numeric(n_t)
  inf   <- numeric(n_t)
  sup   <- numeric(n_t)

  for (t in 1:n_t) {
    vals    <- draws[j, t, ]          # Ns draws del factor j en tiempo t
    media[t] <- mean(vals)
    inf[t]   <- quantile(vals, 0.5 * alpha_sig)        # cuantil 2.5%
    sup[t]   <- quantile(vals, 1 - 0.5 * alpha_sig)    # cuantil 97.5%
  }

  list(media = media, inf = inf, sup = sup)
}

# ============================================================
# 3. Estudio Monte Carlo: cobertura promediada sobre 100 individuos
#
# Para obtener un número estable y
# comparable al 95% nominal, se repite todo el experimento
# (simular α,y → estimar con MMP → calcular cobertura) sobre
# n_individuos réplicas independientes, manteniendo fijos los
# parámetros del modelo (Q, Q1, T, D, a, Z), y se promedia la
# cobertura obtenida en cada una.
# ============================================================
alpha_sig <- 0.05      # significancia
n_individuos <- 100    # num de individuos
Ns           <- 1000   # draws MMP por individuo

cat(sprintf(
  "\n=== Estudio Monte Carlo: cobertura sobre %d individuos (Ns=%d draws c/u) ===\n",
  n_individuos, Ns
))
ptm <- proc.time()

# Matriz de cobertura: una fila por individuo, una columna por factor
cobertura_matriz <- matrix(NA, nrow = n_individuos, ncol = m)

# Se guarda el primer individuo para poder graficarlo de ejemplo más abajo
ejemplo <- NULL

for (i in 1:n_individuos) {

  # --- Simular un nuevo individuo (nueva trayectoria α, y) ---
  sim_i <- modelo_de_simulacion(
    Q     = Q_sim, Q1 = Q1_sim, T_mat = T_sim, D = D_sim,
    a     = a_sim, Z  = Z_sim,  m     = m,     p = p,     n = n
  )

  # --- Estimar sus estados con MMP ---
  resultado_i <- estimar_estados_mmp(
    Y = sim_i$Y, Z = sim_i$Z, T_mat = sim_i$T_mat,
    Q = sim_i$Q, Q1 = sim_i$Q1, D = sim_i$D, a = sim_i$a,
    m = m, n = n, Ns = Ns
  )

  # --- Cobertura del IC 95% de este individuo, por factor ---
  for (j in 1:m) {
    ic         <- calcular_ic(resultado_i$draws, j, alpha_sig)
    alpha_real <- as.vector(sim_i$alpha[j, ])
    cobertura_matriz[i, j] <- mean(alpha_real >= ic$inf & alpha_real <= ic$sup) * 100
  }

  if (i == 1) {
    ejemplo <- list(sim = sim_i, mu_mat = resultado_i$mu_mat, draws = resultado_i$draws)
  }

  if (i %% 10 == 0) cat(sprintf("  Individuo %d/%d completado\n", i, n_individuos))
}

tiempo <- proc.time() - ptm
cat(sprintf("Tiempo total: %.2f segundos\n", tiempo["elapsed"]))

# ============================================================
# 4. Resumen de cobertura promedio (lo que sí debe acercarse a 95%)
# ============================================================
cat(sprintf(
  "\n=== Cobertura del IC al %.0f%%, promediada sobre %d individuos ===\n",
  (1 - alpha_sig) * 100, n_individuos
))
for (j in 1:m) {
  media_j <- mean(cobertura_matriz[, j])
  sd_j    <- sd(cobertura_matriz[, j])
  cat(sprintf(
    "  Factor %d: %.2f%% promedio  (sd entre individuos: %.2f pp, rango: %.1f%% - %.1f%%)\n",
    j, media_j, sd_j, min(cobertura_matriz[, j]), max(cobertura_matriz[, j])
  ))
}
cobertura_global <- mean(cobertura_matriz)
cat(sprintf(
  "\n  Cobertura global (todos los factores, %d individuos): %.2f%%\n",
  n_individuos, cobertura_global
))
cat(sprintf(
  "  (Se espera ~%.0f%% para un IC bien calibrado)\n",
  (1 - alpha_sig) * 100
))

# ============================================================
# 5. Gráfico de ejemplo: un individuo con su banda de IC al 95%
#
# Estructura adaptada de graphstates() en Graphcyd.R:
# - Región sombreada = IC al 95% (polygon)
# - Línea azul punteada = media empírica de los draws
# - Línea negra sólida = estado real (alpha verdadero)
# - Línea roja punteada = media suavizada E[α|y]
#
# Es solo un individuo ilustrativo — la métrica que importa es la
# cobertura promedio de la sección 5, no la de este gráfico aislado.
# ============================================================
tiempo_t <- 1:n

for (j in 1:m) {

  ic <- calcular_ic(ejemplo$draws, j, alpha_sig)

  alpha_real <- as.vector(ejemplo$sim$alpha[j, ])
  dentro     <- mean(alpha_real >= ic$inf & alpha_real <= ic$sup) * 100

  y_min <- min(ic$inf, alpha_real, ejemplo$mu_mat[j, ])
  y_max <- max(ic$sup, alpha_real, ejemplo$mu_mat[j, ])
  y_pad <- (y_max - y_min) * 0.05

  xx <- c(tiempo_t, rev(tiempo_t))
  yy <- c(ic$inf,   rev(ic$sup))

  plot(
    xx, yy,
    type  = "n",
    ylim  = c(y_min - y_pad, y_max + y_pad),
    xlab  = "Tiempo (t)",
    ylab  = bquote(alpha[.(j)]),
    main  = sprintf(
      "Factor %d — individuo de ejemplo (%.1f%% de puntos reales dentro)",
      j, dentro
    )
  )

  polygon(xx, yy,
          col    = adjustcolor("steelblue", alpha.f = 0.2),
          border = "steelblue",
          lty    = "dashed")

  lines(tiempo_t, ic$media,
        col = "steelblue", lwd = 1.5, lty = "dashed")

  lines(tiempo_t, ejemplo$mu_mat[j, ],
        col = "red", lwd = 1.5, lty = "dotted")

  lines(tiempo_t, alpha_real,
        col = "black", lwd = 1.5, lty = "solid")

  legend(
    "topright",
    legend = c(
      "Estado real (α verdadero)",
      sprintf("Media empírica (%d draws)", Ns),
      "Media suavizada E[α|y]",
      sprintf("IC %.0f%%", (1 - alpha_sig) * 100)
    ),
    col    = c("black", "steelblue", "red", adjustcolor("steelblue", 0.4)),
    lwd    = c(1.5, 1.5, 1.5, 8),
    lty    = c("solid", "dashed", "dotted", "solid"),
    cex    = 0.75,
    bty    = "n"
  )
}
