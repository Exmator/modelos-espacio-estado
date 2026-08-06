### Algoritmo de estimación de estados usando el modelo 2###

# ----------------------------------------------------------------
# PASO 1: Calcular Ω y c (Apéndice A simplificado)
# ----------------------------------------------------------------
compute_Omega_c <- function(Y, Z, T_mat, Q, Q1, D, n) {
  # --- Parametros de entrada ---
  # Y = vector px1 de valores observados
  # Z = matriz de coeficientes de las observaciones
  # T_mat = matriz de coeficientes de los estados
  # Q = matriz de covarianzas del error de los estados t
  # Q1 = matriz de covarianzas del error del estado inicial (t = 1)
  # D = matriz de covarianzas del error de las observaciones
  # n = número de períodos de tiempo
  
  
  # Calculo de las inversas
  Q_inv  <- solve(Q)
  Q1_inv <- solve(Q1)
  D_inv  <- solve(D)
  
  # Calculos auxiliares, se repiten varias veces en las ecuaciones, 
  # por lo que se calculan solo una vez
  ZtDZ <- t(Z) %*% D_inv %*% Z          # Z'D⁻¹Z
  TtQT <- t(T_mat) %*% Q_inv %*% T_mat  # T'Q⁻¹T
  
  # Bloques diagonales de Omega
  Omega_diag      <- vector("list", n)
  Omega_diag[[1]] <- ZtDZ + TtQT + Q1_inv     # Ω₁₁
  for (t in 2:(n - 1))
    Omega_diag[[t]] <- ZtDZ + TtQT + Q_inv    # Ωtt, t=2..n-1
  Omega_diag[[n]]   <- ZtDZ + Q_inv           # Ωnn
  
  # Bloque fuera de diagonal: Ω_{t,t+1} = -T'Q⁻¹ (constante en t)
  Omega_off <- -t(T_mat) %*% Q_inv
  
  # Co-vectores: ct = Z'yt para todo t
  c_list <- lapply(1:n, function(t) t(Z) %*% D_inv %*% Y[, t])
  
  list(Omega_diag = Omega_diag, Omega_off = Omega_off, c_list = c_list)
}

# ----------------------------------------------------------------
# PASO 2a: Pre-computación MMP (Result 2.1, forward pass)
# ----------------------------------------------------------------
# Calcula:
#   - U[[t]]: factor de Cholesky de Σt⁻¹ (upper triangular, U't U = Σt⁻¹)
#   - LiO[[t]]: Λt⁻¹ · Ω_{t,t+1}  (para el draw)
#   - m_list[[t]]: medias condicionales mt
# ----------------------------------------------------------------
mmp_precompute <- function(Omega_diag, Omega_off, c_list, n) {
  # --- Parametros de entrada ---
  # Omega_diag = lista de valores de la diagonal de omega
  # Omega_off = valor fuera de la diagonal de omega
  # c_list = vector de co-vectores
  # n = número de períodos de tiempo
  
  U_list   <- vector("list", n)       # Lista de Λt (Cholesky de Σt⁻¹)
  LiO_list <- vector("list", n - 1)   # Lista con valores de [Λt⁻¹ · Ω_{t,t+1}]
  m_list   <- vector("list", n)       # m_t
  
  Sinv <- Omega_diag[[1]]             # Σ₁⁻¹ = Ω₁₁
  
  # Se realiza un ciclo para t = 1,2,...,n
  for (t in 1:n) {
    # Calculo Cholesky para encontrar Λ_t^⊤, donde U = Λ_t^⊤ 
    # tal que t(U) %*% U = Σt⁻¹
    U_list[[t]] <- chol(Sinv)
    
    # Si t no es igual o no es mayor a n
    if (t < n) {
      
      # Computo Λ_t^(-1) · Ω_{t,t+1} usando triangular back-sustitution 
      # donde se resuelve Λ_t · x = Ω_{t,t+1}, 
      # que en el codigo seria t(U) · x = Omega_off
      LiO_list[[t]] <- backsolve(U_list[[t]], Omega_off, transpose = TRUE)
      
      # Encuentro Σ_{t+1}⁻¹ = Ω_{t+1,t+1} - [Ω_{t,t+1}^T · Σ_t · Ω_{t,t+1}]
      # donde Ω_{t,t+1}^T · Σ_t · Ω_{t,t+1} = [Λt⁻¹Ω_{t,t+1}]' [Λt⁻¹Ω_{t,t+1}]
      # que en el codigo seria t(LiO_list[[t]]) %*% LiO_list[[t]] = crossprod(LiO_list[[t]])
      Sinv <- Omega_diag[[t + 1]] - crossprod(LiO_list[[t]])
    }
    
    # Calculo 
    # m1 = (Λt⊤)⁻¹ · Λt⁻¹ · c1 y
    # mt = (Λt⊤)⁻¹ · Λt⁻¹ · (ct - Ω_{t-1,t}' · m_{t-1})
    
    # Defino 
    # si t = 1:     rhs = c1
    # si t != 1:    rhs = ct - Ω_{t-1,t}' · m_{t-1}
    rhs <- if (t == 1) {
      c_list[[1]]
    } else {
      c_list[[t]] - t(Omega_off) %*% m_list[[t - 1]]
    }
    # Resuelvo tmp = Λt⁻¹ · rhs usando triangular back-sustitution (Λt · tmp = rhs)
    tmp         <- backsolve(U_list[[t]], rhs, transpose = TRUE)
    # Resuelvo m_t = (Λt⊤)⁻¹ · tmp usando triangular back-sustitution (Λt⊤ · m_t = tmp)
    m_list[[t]] <- backsolve(U_list[[t]], tmp)                    
  }
  
  # Backward pass para E[α|y] = µ
  mu_list      <- vector("list", n)
  
  # µn = mn
  mu_list[[n]] <- m_list[[n]]
  for (t in (n - 1):1) {
    # µt = mt - Σt · Ω_{t,t+1} · µ_{t+1}
    # Σt · Ω_{t,t+1} = (Λt⊤)⁻¹ · [Λt⁻¹ · Ω_{t,t+1}] = backsolve(U, LiO)
    Sigma_t_Omega <- backsolve(U_list[[t]], LiO_list[[t]])
    mu_list[[t]]  <- m_list[[t]] - Sigma_t_Omega %*% mu_list[[t + 1]]
  }
  
  list(U_list = U_list, LiO_list = LiO_list, m_list = m_list, mu_list = mu_list)
}

# ----------------------------------------------------------------
# PASO 2b: Draw de α | y (backward pass, un draw por llamada)
# ----------------------------------------------------------------
mmp_draw <- function(U_list, LiO_list, m_list, n, m) {
  # --- Parametros de entrada ---
  # U_list = Lista de Λt (Cholesky de Σt⁻¹)
  # LiO_list = Lista con valores de [Λt⁻¹ · Ω_{t,t+1}]
  # m_list = Lista que representa a los m_t
  # n = número de períodos de tiempo
  # m = cantidad de factores del estado
  
  alpha    <- vector("list", n)     # Vector de estados estimado
  
  # t = n, n-1, ..., 1
  for (t in n:1) {
    eps <- rnorm(m)   # εt ~ N(0, Im)
    
    if (t == n) {
      # Para t = n 
      # α_n = m_n + (Λ_n^⊤)⁻¹ · ε_n
      alpha[[n]] <- m_list[[n]] + backsolve(U_list[[n]], eps)
    } else {
      # Para t < n
      # α_t = m_t + (Λ_t^⊤)⁻¹ · (ε_t − [Λ_t⁻¹Ω_{t,t+1}] · α_{t+1}) 
      adj        <- eps - LiO_list[[t]] %*% alpha[[t + 1]]
      alpha[[t]] <- m_list[[t]] + backsolve(U_list[[t]], adj)
    }
  }
  
  do.call(cbind, alpha)  # devuelve matriz m × n
}
