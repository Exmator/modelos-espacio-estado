

# Algoritmo de estimación de estados
mcstates <- function(y,x0,lambda,theta,R,Q,n){
  # ----------------------------------
  #   Objetivo: Algoritmo que retorna el vector nx1 de estados estimados
  # 
  #   Entradas:
  #     y = 
  #     x0 =
  #     lambda =
  #     theta =
  #     R =
  #     Q =
  #     n =
  # 
  #   Salidas:
  #     state = vector nx1 donde cada state_i representa el estado estimado en el tiempo i
  # ----------------------------------
  
  # Definicion de Omega Ω como matriz nxn
  Omega=matrix(0,ncol=n,nrow=n)
  # Definicion de c como vector nx1
  c=array(0,n)
  
  
  #### Calculo de Omega y c ####
  
  # c_1 =  λ_t*R^(−1)*y_1 + Q^(−1)*F_1*x_0
  c[1]=(lambda[1])*(1/R)*y[1]+theta*x0*(Q)^(-1)
  # c_n =  λ_n*R^(−1)*y_n
  c[n]=y[n]*(lambda[n])*(R)^(-1)
  # Ω_(11)= λ_t*R^(−1) + F^T_1*Q^(−1)*F_1 + Q^(−1) (?)
  Omega[1,1]=(lambda[1])*(1/R)+(1/Q)*theta^2+(Q)^(-1)
  # Ω_(21) = −F^T_2*Q^(−1)
  Omega[2,1]=Omega[1,2]=-theta*(1/Q)
  # Ω_(nn)= λ_n*R^(−1) + Q^(−1) + F^T_1*Q^(−1)*F_1 (?)
  Omega[n,n]=(lambda[n])*(1/R)+(1/Q)
  for( t in 2:(n-1)){
    c[t]=(lambda[t])*(1/R)*y[t]
    Omega[t,t]=(lambda[t])*(1/R)+(1/Q)*(theta^2+1)
    Omega[t,t+1]=Omega[t+1,t]=-theta*(1/Q)
  }
  
  #### Calculo de Sigma y mm ####
  
  # Definicion de Sigma Σ como vector nx1
  Sigma = array(0,n)
  # Definicion de Sigma mm como vector nx1
  mm = array(0,n)
  
  
  Sigma[1] = 1/Omega[1,1]
  mm[1]= Sigma[1]*c[1]
  for(t in 2:n){
    Sigma[t] = (1/(Omega[t,t]-Omega[t-1,t]*Sigma[t-1]*Omega[t-1,t]))
    mm[t] = Sigma[t]*(c[t]-Omega[t-1,t]*mm[t-1])
  }
  
  #### Calculo de state ####
  
  # Definicion de state α como vector nx1
  state = array(0,n)
  
  # Definicion de sumastate (?)
  sumastate= array(0,n)
  
  # Dibujo de e_t ∼ N(0,I_m)
  e = rnorm(1)
  # α_n = m_n + (Λ_n^T)^(−1)*e_n 
  state[n]=mm[n]+sqrt(Sigma[n])*e
  # ?
  sumastate[n]=sumastate[n]+state[n]
  
  for(t in (n-1):1){
    e=rnorm(1) 
    state[t]=mm[t]-Sigma[t]*Omega[t,t+1]*state[t+1]+sqrt(Sigma[t])*e
  }
  return(state)
}