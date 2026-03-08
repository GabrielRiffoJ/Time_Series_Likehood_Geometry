
1+0.5^2+0.3^2
0.5-0.5*0.3
-0.3


muestras=c()
for (i  in c(1:500)){
  arma_sim <- arima.sim(n = 10, model = list(ar = c(0.5, -0.3)))
  muestras=rbind(muestras,arma_sim)
  }



resultado=c()

for (i in c(1:500) ){
  Y=muestras[i,]
  n=length(Y)


  fobj<-function(par){
    a0=par[1]
    a1=par[2]
    a2=par[3]
    M=matrix(0,n,n)
    diag(M)=diag(M)+a0^2+a1^2+a2^2
    M[row(M) == col(M) - 1] <- a0*a1+a1*a2
    M[row(M) == col(M) - 2] <- a0*a2
    M[row(M)-1 == col(M)] <- a0*a1+a1*a2
    M[row(M)-2 == col(M)] <- a0*a2
    Minv=chol2inv(chol(M))
    detM=det(M)
    lk=1/2*log(detM)+1/2*t(Y)%*%Minv %*% Y
    return(lk)
  }

  re=optim(c(sd(Y),0,0),fobj)


  a0=re$par[1]
  a1=re$par[2]
  a2=re$par[3]

  resultado=rbind(resultado,c(a0,a1,a2))
}

resultado[1,]

colMeans(resultado)
colMeans(resultado)-c(1,0.5,-0.3)

sd(resultado[,1])
sd(resultado[,2])
sd(resultado[,3])
mean(resultado[,1]-1)
mean(resultado[,2]-0.5)
mean(resultado[,3]+0.3)

sum((resultado[9,]-c(1,0.5,-0.3))^2)

sum((c(0.1124229,0.0395593,-0.45528))^2)






#save samples for HomotopyContinuation

write.table(muestras, file = "samples.txt", row.names = FALSE, col.names = TRUE, sep = "\t")




#result=c()

a0
a1
a2

sum(abs(c(a0-1,a1-0.5, a2+0.3)))

a0^2+a1^2+a2^2
a0*a1+a1*a2
a0*a2





