
set.seed(1234)
muestras=c()
for (i in c(1:500)){
  arma_sim <- arima.sim(n = 15, model = list(ar = 0.5))
  muestras=rbind(muestras,arma_sim)
}

resultado=c()

for (i in c(1:500) ){
  Y=muestras[i,]
  n=length(Y)

  fobj<-function(par){
    phi=par[1]
    s2=par[2]
    if (abs(phi)>=0.999 || s2<=0) return(1e10)
    M=matrix(0,n,n)
    M=s2/(1-phi^2)*phi^abs(row(M)-col(M))
    Minv=chol2inv(chol(M))
    detM=det(M)
    lk=1/2*log(detM)+1/2*t(Y)%*%Minv %*% Y
    return(lk)
  }

  re=optim(c(0,var(Y)),fobj)

  phi=re$par[1]
  s2=re$par[2]

  resultado=rbind(resultado,c(phi,s2))
}

resultado[1,]

colMeans(resultado)
colMeans(resultado)-c(0.5,1)

sd(resultado[,1])
sd(resultado[,2])
mean(resultado[,1]-0.5)
mean(resultado[,2]-1)

write.table(muestras, file = "samples.txt", row.names = FALSE, col.names = TRUE, sep = "\t")





