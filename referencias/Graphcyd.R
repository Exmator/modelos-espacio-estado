#------- Graficos

#-------Varianzas
graphvar=function(index,sigma,Q,sigmav,Qv){
  M = matrix(c(1,2,3,3), byrow=T, ncol=2)
  layout(M)
  plot(Q[index,1],lty=1,type="l", ylab=expression(Q[11]))
  abline(h=Qv[1,1],col="red")
  plot(Q[index,4],lty=1,type="l",ylab=expression(Q[22]))
  abline(h=Qv[2,2],col="red")
  plot(sigma[index],lty=1,type="l",ylab=expression(sigma^2))
  abline(h=sigmav,col="red") 
}

graphD=function(index,D,Dv){
  M = matrix(c(1,2,3,3), byrow=T, ncol=2)
  layout(M)
  plot(D[index,1],lty=1,type="l",ylab=expression(D[11]))
  abline(h=Dv[1,1],col="red")
  plot(D[index,2],lty=1,type="l",ylab=expression(D[22]))
  abline(h=Dv[2,2],col="red")
  plot(D[index,3],lty=1,type="l",ylab=expression(D[55]))
  abline(h=Dv[5,5],col="red")
}

#---- Parametros
graphpop=function(index,thetapop,thetav){
  M = matrix(c(1,2,3,3), byrow=T, ncol=2)
  layout(M)
  plot(thetapop[index,1],type="l",ylab=expression(lambda))
  abline(h=thetav[1],col="red")
  plot(thetapop[index,2],type="l",ylab=expression(alpha))
  abline(h=thetav[2],col="red")
  plot(thetapop[index,3],type="l",ylab=expression(N))
  abline(h=thetav[5],col="red")
}

graphind=function(index,thetas1,thetav){
  M = matrix(c(1,2,3,3), byrow=T, ncol=2)
  layout(M)
  plot(thetas1[index,1],type="l",ylab=expression(lambda[i]))
  abline(h=thetav[1],col="red")
  plot(thetas1[index,2],type="l",ylab=expression(alpha[i]))
  abline(h=thetav[2],col="red")
  plot(thetas1[index,3],type="l",ylab=expression(N[i]))
  abline(h=thetav[5],col="red")
}

#------ States
graphstates=function(index,tem,xsam1,xsam2,gsam,dgsam,x,pc){
  desvio1=0
  desvio2=0
  desvio3=0
  desvio4=0
  
  media1=0
  media2=0
  media3=0
  media4=0
  
  inf1=0
  inf2=0
  inf3=0
  inf4=0
  
  sup1=0
  sup2=0
  sup3=0
  sup4=0
  
  
  for(i in 1:n){desvio1[i]=sd(xsam1[i,index])
  media1[i]=mean(xsam1[i,index])
  inf1[i]=quantile(xsam1[i,index],0.025)
  sup1[i]=quantile(xsam1[i,index],0.975)}
  for(i in 1:n){desvio2[i]=sd(xsam2[i,index])
  media2[i]=mean(xsam2[i,index])
  inf2[i]=quantile(xsam2[i,index],0.025)
  sup2[i]=quantile(xsam2[i,index],0.975)}
  for(i in 1:n){desvio3[i]=sd(gsam[i,index])
  media3[i]=mean(gsam[i,index])
  inf3[i]=quantile(gsam[i,index],0.025)
  sup3[i]=quantile(gsam[i,index],0.975)}
  for(i in 1:n){desvio4[i]=sd(dgsam[i,index])
  media4[i]=mean(dgsam[i,index])
  inf4[i]=quantile(dgsam[i,index],0.025)
  sup4[i]=quantile(dgsam[i,index],0.975)}
  
  
  ICi1=media1-2.5*desvio1
  ICs1=media1+2.5*desvio1
  
  ICi2=media2-2.5*desvio2
  ICs2=media2+2.5*desvio2
  
  ICi3=media3-2.5*desvio3
  ICs3=media3+2.5*desvio3
  
  
  
  time=tem      
  par(mfrow=c(2,2))
  xx <- c(time, rev(time))
  yy1 <- c(inf1,rev(sup1))
  plot1=plot(xx, yy1,main="T+T*", type = "n", xlab = "Day", ylab = "T-cell")
  polygon(xx, yy1, col = "white",lty = "longdash", border = "blue",xlab="Time",ylab="ht")
  lines(time,media1, lty = "longdash",col="blue")
  lines(time,x[(4*pc-3),], lty = "solid",col="black")
  
  yy2 <- c(inf2,rev(sup2))
  plot1=plot(xx, yy2,main="V", type = "n", xlab = "Day", ylab = "Viral Load ")
  polygon(xx, yy2, col = "white", lty = "longdash",border = "blue",xlab="Time",ylab="ht")
  lines(time,media2, lty = "longdash",col="blue")
  lines(time,x[(4*pc-2),], lty = "solid",col="black")
  
  
  yy3 <- c(inf3,rev(sup3))
  plot1=plot(xx, yy3,main="g(t)", type = "n", xlab = "t", ylab = "g(t) ")
  polygon(xx, yy3, col = "white",lty = "longdash", border = "blue",xlab="Time",ylab="ht")
  lines(time,media3, lty = "longdash",col="blue")
  lines(time,x[(4*pc-1),], lty = "solid",col="black")
  
  yy4 <- c(inf4,rev(sup4))
  plot1=plot(xx, yy4,main="dg(t)/dt", type = "n", xlab = "t", ylab = "g(t) ")
  polygon(xx, yy4, col = "white",lty = "longdash", border = "blue",xlab="Time",ylab="ht")
  lines(time,media4, lty = "longdash",col="blue")
  lines(time,x[(4*pc),], lty = "solid",col="black")
  
}

ICstates=function(index,xsam1,xsam2,gsam,alpha){
  n=length(xsam1[,1])
  desvio1=0
  desvio2=0
  desvio3=0
  
  media1=0
  media2=0
  media3=0
  
  inf1=0
  inf2=0
  inf3=0
  
  sup1=0
  sup2=0
  sup3=0
  
  
  for(i in 1:n){desvio1[i]=sd(xsam1[i,index])
  media1[i]=mean(xsam1[i,index])
  inf1[i]=quantile(xsam1[i,index],0.5*alpha)
  sup1[i]=quantile(xsam1[i,index],1-0.5*alpha)}
  for(i in 1:n){desvio2[i]=sd(xsam2[i,index])
  media2[i]=mean(xsam2[i,index])
  inf2[i]=quantile(xsam2[i,index],0.5*alpha)
  sup2[i]=quantile(xsam2[i,index],1-0.5*alpha)}
  for(i in 1:n){desvio3[i]=sd(gsam[i,index])
  media3[i]=mean(gsam[i,index])
  inf3[i]=quantile(gsam[i,index],0.5*alpha)
  sup3[i]=quantile(gsam[i,index],1-0.5*alpha)}
  
  
  outi1=which(x[(4*pc-3),]<inf1)
  outs1=which(x[(4*pc-3),]>sup1)
  
  outi2=which(x[(4*pc-2),]<inf2)
  outs2=which(x[(4*pc-2),]>sup2)
  
  outi3=which(x[(4*pc-1),]<inf3)
  outs3=which(x[(4*pc-1),]>sup3)
  
  return(list(IC1=c(outi1,outs1),IC2=c(outi2,outs2),IC3=c(outi3,outs3)))
}


