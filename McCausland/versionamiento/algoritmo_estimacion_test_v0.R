# ================================================================
# EJECUCIÓN Y VERIFICACIÓN DEL ALGORITMO DE ESTIMACION
# ================================================================
source("algoritmo_estimacion.R")
source("algoritmo_simulacion.R")

sim <- modelo_de_simulacion(m = 4, p = 10, n = 1000)

# Calcular Ω y c
omega_c <- compute_Omega_c(
  Y = sim$Y, Z = sim$Z, T_mat = sim$T,
  Q = sim$Q, Q1 = sim$Q1, D = sim$D,
  n = sim$n
)

# Pre-computación
precomp <- mmp_precompute(
  Omega_diag = omega_c$Omega_diag,
  Omega_off  = omega_c$Omega_off,
  c_list     = omega_c$c_list,
  n = sim$n
)

# Un draw de α | y
alpha_draw <- mmp_draw(
  U_list = precomp$U_list, 
  LiO_list = precomp$LiO_list, 
  m_list = precomp$m_list, 
  n = sim$n, 
  m = sim$m
  )

# Media suavizada E[α | y]
mu_mat <- matrix(unlist(precomp$mu_list), nrow = sim$m)  # m × n

# ---- Verificación: comparar media suavizada con estados verdaderos ----
cat("Correlación entre E[α|y] y α verdadero por factor:\n")

for (j in 1:sim$m) {
  cat(sprintf("  Factor %d: %.4f\n", j, cor(mu_mat[j, ], sim$alpha[j, ])))
}

cat("\nDraw único — primeras 5 columnas (factor 1):\n")
print(round(alpha_draw[1, 1:5], 4))

cat("\nMedia suavizada — primeras 5 columnas (factor 1):\n")
print(round(mu_mat[1, 1:5], 4))

# Múltiples draws (para importance sampling, por ejemplo)
Ns     <- 50
draws  <- array(0, dim = c(sim$m, sim$n, Ns))
for (s in 1:Ns) draws[, , s] <- mmp_draw(U_list = precomp$U_list, 
                                         LiO_list = precomp$LiO_list, 
                                         m_list = precomp$m_list, 
                                         n = sim$n, 
                                         m = sim$m
                                         )

cat("\nMedia empírica de", Ns, "draws (factor 1, t=1:5) vs media suavizada:\n")
empirica <- apply(draws[1, 1:5, ], 1, mean)
cat("  Empírica:   ", round(empirica, 4), "\n")
cat("  Suavizada:  ", round(mu_mat[1, 1:5], 4), "\n")