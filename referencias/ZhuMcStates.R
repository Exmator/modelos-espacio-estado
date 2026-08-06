library(MASS)
library(Matrix)
library(LearnBayes)
library(mvtnorm)

FF=function(theta,h){ #theta=(si,ai,bi,ri,Ni)
  f1= (1-h*theta[3])
  f2=theta[5]*theta[3]*h
  f3=(1-theta[4]*h)
  FF=matrix(c(f1,0,f2,f3),ncol=2,byrow=T)
  return(FF)
}

cc=function(theta,h){
  rho=matrix(c(theta[1]*h,0),nrow=2,ncol=1,byrow=T)
  return(rho)
}

qq=function(theta,h){
  rho=matrix(c(h*(theta[3]-theta[2]),-theta[5]*theta[3]*h),nrow=2,ncol=1,byrow=T)
  return(rho)
}

gt=function(t,h){
  gt=50+ 50*cos(0.8*pi*t*h)
  return(gt)}

dgt=function(t,h){
  dgt=-40*pi*sin(0.8*pi*t*h)
  return(dgt)}

#------ simulation

sim.messmtvp=function(x0,Z,theta,sigmav,Q,D,m,n,r,h){
  x=matrix(0,nrow=2*m,ncol=n)
  y=matrix(0,nrow=n,ncol=m)
  thetas=matrix(0,nrow=m,ncol=r)
  for(i in 1:m){
    thetas[i,]=theta+sqrt(diag(D))*rnorm(r)
    x[(2*i-1):(2*i),1]=cc(thetas[i,],h)+FF(thetas[i,],h)%*%x0[i,]+gt(0,h)*qq(thetas[i,],h)+mvrnorm(1,mu=rep(0,2),Sigma=Q)
    y[1,i]=(Z%*%x[(2*i-1):(2*i),1])+ sqrt(sigmav)*rnorm(1)
    for( t in 2:n){
      x[(2*i-1):(2*i),t]=cc(thetas[i,],h)+FF(thetas[i,],h)%*% x[(2*i-1):(2*i),t-1]+gt(t-1,h)*qq(thetas[i,],h)+mvrnorm(1,mu=rep(0,2),Sigma=Q)
      y[t,i]=(Z%*%x[(2*i-1):(2*i),t])+ sqrt(sigmav)*rnorm(1)
    }
  }
  return(list(x=x,y=y,thetas=thetas))
}

#---------- Funciones
FFtil=function(theta,h){ #theta=(si,ai,bi,ri,Ni)
  f1= (1-h*theta[3])
  f2=h*(theta[3]-theta[2])
  f3=theta[3]*theta[5]*h
  f4=(1-theta[4]*h)
  f5=-h*theta[5]*theta[3]
  FF=matrix(c(f1,0,f2,0,f3,f4,f5,0,0,0,1,h,0,0,0,1),ncol=4,byrow=T)
  return(FF)
}
cc=function(theta,h){
  rho=matrix(c(theta[1]*h,0),nrow=2,ncol=1,byrow=T)
  return(rho)
}

R=function(h){
  #r=c(h,0.5*h^2,0.5*h^2,(h^3)/3) #Como en Wecker
  r=c((h^3)/3,0.5*h^2,0.5*h^2,h)  #Como en Zhu 
  R=matrix(r,ncol=2,byrow=T) 
  return(R)
}

#------ Mc States
State=Data$x
x=matrix(0,nrow=4*m,ncol=n)
for(i in 1:m){
  x[(4*i-3):(4*i),]=rbind(State[(2*i-1),],State[(2*i),],gt(1:n,h),dgt(1:n,h))
}
i=1
y=Data$y[,i]
thetas=Data$thetas[i,]
Q=Qv
g0=c(100,0)
x0=x0v[i,]
RR=R(h)
lambda=3.9*10^(-5)
tem=h*(1:n)
mcstatestvp=function(y,Z,x0,g0,thetas,Q,sigmav,RR,n,h,lambda){
  x00=matrix(c(x0,g0),ncol=1,byrow=T)
  ZZ=matrix(c(Z,0,0),ncol=4,byrow=T)
  Wb=matrix(c(cc(thetas,h),0,0),ncol=1,byrow=T)
  TT=FFtil(thetas,h)
  a1=Wb+TT%*%(x00) 
  P1=as.matrix(bdiag(Q,RR/lambda))
  A11=1/sigmav
  A22=solve(as.matrix(bdiag(Q,RR/lambda)))
  Omega11=t(ZZ)%*%(ZZ)*A11+t(TT)%*%A22%*%(TT)+solve(P1)
  Omegat1=-t(TT)%*%A22
  Omegatt=t(ZZ)%*%(ZZ)*A11+t(TT)%*%A22%*%(TT)+A22
  c=matrix(0,ncol=4,nrow=n)
  c[1,]=as.vector(t(ZZ)*A11*y[1]-t(TT)%*%A22%*%Wb+solve(P1)%*%a1)
  for(t in 2:(n-1)){
    c[t,]=as.vector(t(ZZ)*A11*y[t]-t(TT)%*%A22%*%Wb+A22%*%Wb)
  }
  c[n,]=as.vector(t(ZZ)*A11*y[n]+A22%*%Wb)
  Omegann=A11*t(ZZ)%*%ZZ+A22 
  
  Sigma=matrix(0,nrow=(4*n),ncol=4)
  mm=matrix(0,ncol=4,nrow=n)
  state = matrix(0,nrow=4,ncol=n)
  
  Sigma[1:4,]=as.matrix(solve(Omega11))
  mm[1,]=Sigma[1:4,]%*%c[1,]
  
  for(t in 2:(n-1)){
    Sigma[(4*t-3):(4*t),]=as.matrix(solve(Omegatt-t(Omegat1)%*%Sigma[(4*t-7):(4*t-4),]%*%(Omegat1)))
    mm[t,]=as.vector(Sigma[(4*t-3):(4*t),]%*%(c[t,]-t(Omegat1)%*%mm[t-1,]))
  }
  Sigma[(4*n-3):(4*n),]=as.matrix(solve(Omegann-t(Omegat1)%*%Sigma[(4*n-7):(4*n-4),]%*%(Omegat1)))
  mm[n,]=as.vector(Sigma[(4*n-3):(4*n),]%*%(c[n,]-t(Omegat1)%*%mm[n-1,]))
  
  state[,n]=mvrnorm(1,mm[n,],Sigma=Sigma[(4*n-3):(4*n),])
  for(t in (n-1):1){
    mu=mm[t,]-as.vector(Sigma[(4*t-3):(4*t),]%*%Omegat1%*%state[,(t+1)])
    state[,t]=mvrnorm(1,mu,Sigma=Sigma[(4*t-3):(4*t),])
  }
  #par(mfrow=c(2,2))
  #plot(tem,state[1,],type="l")
  #lines(tem,x[(4*i-3),],col="red")
  #plot(tem,state[2,],type="l")
  #lines(tem,x[(4*i-2),],col="red")
  #plot(tem,state[3,],type="l")
  #lines(tem,x[(4*i-1),],col="red")
  #plot(tem,state[4,],type="l")
  #lines((1:n)*h,x[(4*i),],col="red")
  
  return(state)}

#----- Simulated Data
h=0.06
n=50
m=2
Z=matrix(c(0,1),ncol=2,byrow=T)
sigmav=2*10^6
Qv=diag(c(0.4^2,5*10^3))
Dv=diag(c(5^2,0.01^2,0.1^2,0.5^2,100))
thetav=c(30,0.1,0.5,3,1000)
x0v=mvrnorm(m,mu=c(150,5*10^4),Sigma=0.0025*Qv) 
Data=sim.messmtvp(x0v,Z,thetav,sigmav,Qv,Dv,m,n,length(thetav),h)
par(mfrow=c(1,2))
plot((1:n)*h,Data$y[,1],type="l")
plot((1:n)*h,Data$y[,2],type="l")

par(mfrow=c(2,2))
plot((1:n)*h,Data$x[1,],type="p",pch = 18)
plot((1:n)*h,Data$x[2,],type="p",pch = 18)
plot((1:n)*h,Data$y[,1],type="p",pch = 18)
plot((1:n)*h,gt(1:n,h),type="l")


#------ objetos
sam=5000
xsam11=matrix(0,nrow=n,ncol=sam)
xsam12=matrix(0,nrow=n,ncol=sam)
gsam1=matrix(0,nrow=n,ncol=sam)
dgsam1=matrix(0,nrow=n,ncol=sam)
xsam21=matrix(0,nrow=n,ncol=sam)
xsam22=matrix(0,nrow=n,ncol=sam)
gsam2=matrix(0,nrow=n,ncol=sam)
dgsam2=matrix(0,nrow=n,ncol=sam)

#----Valores Iniciales
State=Data$x
x=matrix(0,nrow=4*m,ncol=n)
for(i in 1:m){
  x[(4*i-3):(4*i),]=rbind(State[(2*i-1),],State[(2*i),],gt(1:n,h),dgt(1:n,h))
}
xupdn=0.95*x
x0=x0v
g0=matrix(c(100,0),nrow=m,ncol=2,byrow=T)
Qsim=Qv
thetasi=Data$thetas
yobs=Data$y
ptm <- proc.time()
for(iter in 1:sam){
  for(i in 1:m){
    xupdn[(4*i-3):(4*i),]=mcstatestvp(yobs[,i],Z,x0[i,],g0[i,],thetasi[i,],Qsim,sigmav,R(h),n,h,3.9*10^(-5))
    }
  print(iter)
  xsam11[,iter]=xupdn[1,]
  xsam12[,iter]=xupdn[2,]
  gsam1[,iter]=xupdn[3,]
  dgsam1[,iter]=xupdn[4,]
  xsam21[,iter]=xupdn[5,]
  xsam22[,iter]=xupdn[6,]
  gsam2[,iter]=xupdn[7,]
  dgsam2[,iter]=xupdn[8,]
}
tiempo=proc.time()-ptm

index=seq(sam/2,sam,1)
tem=h*(1:n)
xsam1=xsam21
xsam2=xsam22
gsam=gsam2
dgsam=dgsam2
pc=2
graphstates(index,tem,xsam1,xsam2,gsam,dgsam,x,pc)
#------ IC States
ICstates(index,xsam1,xsam2,gsam,0.01)
