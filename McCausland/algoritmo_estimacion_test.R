# ============================================================
# Test del algoritmo de estimación MMP
# McCausland et al. (2011) - Simulation smoothing for state-space models
#
# Compara estados reales vs estimados con IC al 95% (alpha=0.05)
# usando 1000 draws del algoritmo MMP.
# Caso simple: m=2, p=2 para reducir tiempo de cómputo.
# ============================================================

source("algoritmo_simulacion.R")
source("algoritmo_estimacion.R")

set.seed(42)

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
# 2. Simular datos
# ============================================================
cat("\n=== Simulando datos ===\n")
sim <- modelo_de_simulacion(
  Q     = Q_sim,
  Q1    = Q1_sim,
  T_mat = T_sim,
  D     = D_sim,
  a     = a_sim,
  Z     = Z_sim,
  m     = m,
  p     = p,
  n     = n
)
cat("Datos simulados correctamente.\n")
cat("Dimensión Y (observaciones):", nrow(sim$Y), "×", ncol(sim$Y), "\n")
cat("Dimensión alpha (estados):  ", nrow(sim$alpha), "×", ncol(sim$alpha), "\n")

# ============================================================
# 3. Estimar estados con el algoritmo MMP (1000 draws)
# ============================================================
Ns <- 1000
cat(sprintf("\n=== Estimando estados con MMP (%d draws) ===\n", Ns))
ptm <- proc.time()

resultado <- estimar_estados_mmp(
  Y     = sim$Y,
  Z     = sim$Z,
  T_mat = sim$T_mat,
  Q     = sim$Q,
  Q1    = sim$Q1,
  D     = sim$D,
  a     = sim$a,
  m     = m,
  n     = n,
  Ns    = Ns
)

tiempo <- proc.time() - ptm
cat(sprintf("Tiempo de ejecución: %.2f segundos\n", tiempo["elapsed"]))

mu_mat <- resultado$mu_mat   # m × n  — media suavizada
draws  <- resultado$draws    # m × n × Ns — draws

# ============================================================
# 4. Verificación numérica
# ============================================================
cat("\n=== Verificación: correlación E[α|y] vs α real ===\n")
for (j in 1:m) {
  r <- cor(as.vector(mu_mat[j, ]), as.vector(sim$alpha[j, ]))
  cat(sprintf("  Factor %d: %.4f\n", j, r))
}

# ============================================================
# 5. Construir intervalos de confianza al 95% (alpha = 0.05)
#
# Adaptado de Graphcyd.R (ICstates) y ZhuMcStates.R:
# Para cada factor j y cada tiempo t, los Ns draws forman
# una distribución empírica. Se extraen los cuantiles 2.5% y 97.5%
# como límites del IC, y la media empírica como estimación central.
#
# draws[j, t, ] es el vector de Ns valores del factor j en el tiempo t.
# ============================================================
alpha_sig <- 0.05   # significancia del intervalo de confianza

calcular_ic <- function(draws, j, alpha_sig) {
  # Para el factor j, calcula media y límites del IC para cada t
  # draws: array m×n×Ns
  n_t <- dim(draws)[2]
  Ns  <- dim(draws)[3]

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
# 6. Gráficos: estado real vs estimado con IC al 95%
#
# Estructura adaptada de graphstates() en Graphcyd.R:
# - Región sombreada = IC al 95% (polygon)
# - Línea azul punteada = media empírica de los draws
# - Línea negra sólida = estado real (alpha verdadero)
# - Línea roja punteada = media suavizada E[α|y]
# ============================================================
tiempo_t <- 1:n

# Un gráfico por factor
for (j in 1:m) {

  ic <- calcular_ic(draws, j, alpha_sig)

  # Porcentaje de puntos reales dentro del IC
  alpha_real  <- as.vector(sim$alpha[j, ])
  dentro       <- mean(alpha_real >= ic$inf & alpha_real <= ic$sup) * 100

  # Rango y del gráfico
  y_min <- min(ic$inf, alpha_real, mu_mat[j, ])
  y_max <- max(ic$sup, alpha_real, mu_mat[j, ])
  y_pad <- (y_max - y_min) * 0.05

  # Vectores para el polígono del IC (izquierda → derecha → derecha → izquierda)
  xx <- c(tiempo_t, rev(tiempo_t))
  yy <- c(ic$inf,   rev(ic$sup))

  plot(
    xx, yy,
    type  = "n",
    ylim  = c(y_min - y_pad, y_max + y_pad),
    xlab  = "Tiempo (t)",
    ylab  = bquote(alpha[.(j)]),
    main  = sprintf(
      "Factor %d — IC %.0f%% (%.1f%% de puntos reales dentro)",
      j, (1 - alpha_sig) * 100, dentro
    )
  )

  # Región sombreada del IC
  polygon(xx, yy,
          col    = adjustcolor("steelblue", alpha.f = 0.2),
          border = "steelblue",
          lty    = "dashed")

  # Media empírica de los draws
  lines(tiempo_t, ic$media,
        col = "steelblue", lwd = 1.5, lty = "dashed")

  # Media suavizada E[α|y]
  lines(tiempo_t, mu_mat[j, ],
        col = "red", lwd = 1.5, lty = "dotted")

  # Estado real
  lines(tiempo_t, alpha_real,
        col = "black", lwd = 1.5, lty = "solid")

  # Leyenda
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

# ============================================================
# 7. Resumen de cobertura del IC
# ============================================================
cat(sprintf("\n=== Cobertura del IC al %.0f%% ===\n", (1 - alpha_sig) * 100))
for (j in 1:m) {
  ic           <- calcular_ic(draws, j, alpha_sig)
  alpha_real   <- as.vector(sim$alpha[j, ])
  cobertura    <- mean(alpha_real >= ic$inf & alpha_real <= ic$sup) * 100
  cat(sprintf("  Factor %d: %.1f%% de puntos dentro del IC\n", j, cobertura))
}
cat(sprintf(
  "\n  (Se espera ~%.0f%% para un IC bien calibrado)\n",
  (1 - alpha_sig) * 100
))
