# ============================================================
# Algoritmo tomado de:
# McCausland et al. (2011) - Simulation smoothing for state–space models
# Sección 2: Precision-based methods for simulation smoothing
# ============================================================

# ----------------------------------------------------------------
# PASO 1: Calcular Ω y c (Apéndice A simplificado)
# ----------------------------------------------------------------
compute_Omega_c <- function(Y, Z, T_mat, Q, Q1, D, a, n) {
  # --- Parámetros de entrada ---
  # Y     = matriz p×n de observaciones
  # Z     = matriz p×m de factor loadings
  # T_mat = matriz m×m de transición de estados
  # Q     = matriz m×m de covarianza del error de los estados t
  # Q1    = matriz m×m de covarianza del error del estado inicial (t=1)
  # D     = matriz p×p de covarianza del error de las observaciones
  # a     = media del estado inicial α₁ (m×1)
  # n     = número de períodos de tiempo

  # Cálculo de las inversas
  Q_inv  <- solve(Q)
  Q1_inv <- solve(Q1)
  D_inv  <- solve(D)

  # Cálculos auxiliares (constantes en t, se calculan una sola vez)
  ZtDZ <- t(Z) %*% D_inv %*% Z          # Z'D⁻¹Z  (m×m)
  TtQT <- t(T_mat) %*% Q_inv %*% T_mat  # T'Q⁻¹T  (m×m)

  # Bloques diagonales de Omega
  Omega_diag      <- vector("list", n)
  Omega_diag[[1]] <- ZtDZ + TtQT + Q1_inv    # Ω₁₁
  for (t in 2:(n - 1))
    Omega_diag[[t]] <- ZtDZ + TtQT + Q_inv   # Ωtt,  t = 2,...,n-1
  Omega_diag[[n]]   <- ZtDZ + Q_inv          # Ωnn

  # Bloque fuera de diagonal: Ω_{t,t+1} = -T'Q⁻¹ (constante en t)
  Omega_off <- -t(T_mat) %*% Q_inv

  # Co-vectores: c1 = Z'D⁻¹y1 + Q1⁻¹a,  ct = Z'D⁻¹yt para t = 2,...,n
  c_list      <- lapply(1:n, function(t) t(Z) %*% D_inv %*% Y[, t])
  c_list[[1]] <- c_list[[1]] + Q1_inv %*% a

  list(Omega_diag = Omega_diag, Omega_off = Omega_off, c_list = c_list)
}

# ----------------------------------------------------------------
# PASO 2a: Pre-computación MMP (Result 2.1 — forward pass)
# Calcula:
#   U_list[[t]]   : Cholesky de Σt⁻¹  (triangular superior, U = Λt⊤)
#   LiO_list[[t]] : Λt⁻¹ · Ω_{t,t+1}
#   m_list[[t]]   : medias condicionales mt
#   mu_list[[t]]  : E[αt | y]  (media suavizada, backward pass)
# ----------------------------------------------------------------
mmp_precompute <- function(Omega_diag, Omega_off, c_list, n) {
  # --- Parámetros de entrada ---
  # Omega_diag = lista de n bloques diagonales de Omega
  # Omega_off  = bloque fuera de diagonal (constante en t)
  # c_list     = lista de n co-vectores ct
  # n          = número de períodos de tiempo

  U_list   <- vector("list", n)      # Λt⊤ (Cholesky de Σt⁻¹)
  LiO_list <- vector("list", n - 1)  # Λt⁻¹ · Ω_{t,t+1}
  m_list   <- vector("list", n)      # mt

  Sinv <- Omega_diag[[1]]            # Σ₁⁻¹ = Ω₁₁

  # Se realiza un ciclo para t = 1,2,...,n
  for (t in 1:n) {

    # Paso 1: Cholesky de Σt⁻¹ = Λt · Λt⊤
    # chol() devuelve U = Λt⊤ (triangular superior)
    U_list[[t]] <- chol(Sinv)

    if (t < n) {
      # Paso 2: Λt⁻¹ · Ω_{t,t+1}
      # Resuelve el sistema triangular Λt · x = Ω_{t,t+1}
      LiO_list[[t]] <- backsolve(U_list[[t]], Omega_off, transpose = TRUE)

      # Paso 3: Actualizar Σ_{t+1}⁻¹
      # Σ_{t+1}⁻¹ = Ω_{t+1,t+1} - [Λt⁻¹Ω_{t,t+1}]⊤ [Λt⁻¹Ω_{t,t+1}]
      Sinv <- Omega_diag[[t + 1]] - crossprod(LiO_list[[t]])
    }

    # Paso 4: Calcular mt
    # rhs = ct           si t = 1
    # rhs = ct - Ω_{t-1,t}⊤ · m_{t-1}  si t > 1
    rhs <- if (t == 1) {
      c_list[[1]]
    } else {
      c_list[[t]] - t(Omega_off) %*% m_list[[t - 1]]
    }
    # mt = (Λt⊤)⁻¹ · Λt⁻¹ · rhs  (dos sistemas triangulares)
    tmp         <- backsolve(U_list[[t]], rhs, transpose = TRUE)
    m_list[[t]] <- backsolve(U_list[[t]], tmp)
  }

  # --- Backward pass: E[α|y] = µ ---
  # µn = mn,  µt = mt - Σt·Ω_{t,t+1}·µ_{t+1}
  mu_list      <- vector("list", n)
  mu_list[[n]] <- m_list[[n]]
  for (t in (n - 1):1) {
    Sigma_t_Omega <- backsolve(U_list[[t]], LiO_list[[t]])
    mu_list[[t]]  <- m_list[[t]] - Sigma_t_Omega %*% mu_list[[t + 1]]
  }

  list(U_list = U_list, LiO_list = LiO_list, m_list = m_list, mu_list = mu_list)
}

# ----------------------------------------------------------------
# PASO 2b: Draw de α | y  (backward draw, un draw por llamada)
# ----------------------------------------------------------------
mmp_draw <- function(U_list, LiO_list, m_list, n, m) {
  # --- Parámetros de entrada ---
  # U_list   = lista de Λt⊤ (Cholesky de Σt⁻¹)
  # LiO_list = lista de Λt⁻¹ · Ω_{t,t+1}
  # m_list   = lista de mt
  # n        = número de períodos de tiempo
  # m        = dimensión del estado

  alpha <- vector("list", n)

  for (t in n:1) {
    eps <- rnorm(m)   # εt ~ N(0, Im)

    if (t == n) {
      # αn = mn + (Λn⊤)⁻¹ · εn
      alpha[[n]] <- m_list[[n]] + backsolve(U_list[[n]], eps)
    } else {
      # αt = mt + (Λt⊤)⁻¹ · (εt - [Λt⁻¹Ω_{t,t+1}] · α_{t+1})
      adj        <- eps - LiO_list[[t]] %*% alpha[[t + 1]]
      alpha[[t]] <- m_list[[t]] + backsolve(U_list[[t]], adj)
    }
  }

  matrix(unlist(alpha), nrow = m)   # devuelve matriz m × n
}

# ----------------------------------------------------------------
# FUNCIÓN PRINCIPAL: Ejecuta los tres pasos del algoritmo MMP
# ----------------------------------------------------------------
estimar_estados_mmp <- function(Y, Z, T_mat, Q, Q1, D, a, m, n, Ns = 1) {
  # --- Parámetros de entrada ---
  # Y     = matriz p×n de observaciones
  # Z     = matriz p×m de factor loadings
  # T_mat = matriz m×m de transición de estados
  # Q     = matriz m×m de covarianza del error de los estados
  # Q1    = matriz m×m de covarianza del error del estado inicial
  # D     = matriz p×p de covarianza del error de observaciones
  # a     = media del estado inicial α₁ (m×1)
  # m     = dimensión del estado (número de factores)
  # n     = número de períodos de tiempo
  # Ns    = número de draws a realizar (por defecto 1000)
  #
  # --- Retorna ---
  # mu_mat  : matriz m×n con la media suavizada E[α|y]
  # draws   : array m×n×Ns con los Ns draws de α|y

  # Paso 1: Calcular Ω y c
  omega_c <- compute_Omega_c(
    Y = Y, Z = Z, T_mat = T_mat,
    Q = Q, Q1 = Q1, D = D, a = a, n = n
  )

  # Paso 2a: Pre-computación (forward + backward pass)
  precomp <- mmp_precompute(
    Omega_diag = omega_c$Omega_diag,
    Omega_off  = omega_c$Omega_off,
    c_list     = omega_c$c_list,
    n = n
  )

  # Media suavizada E[α|y]
  mu_mat <- matrix(unlist(precomp$mu_list), nrow = m)

  # Paso 2b: Ns draws de α|y
  draws <- array(0, dim = c(m, n, Ns))
  for (s in 1:Ns) {
    draws[, , s] <- mmp_draw(
      U_list   = precomp$U_list,
      LiO_list = precomp$LiO_list,
      m_list   = precomp$m_list,
      n = n, m = m
    )
  }

  list(mu_mat = mu_mat, draws = draws)
}
