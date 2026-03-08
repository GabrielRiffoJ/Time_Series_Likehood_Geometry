0.5^2+0.3^2+1
muestras2=c()
for (i  in c(1:500)){
  arma_sim <- arima.sim(n = 500, model = list(ma = c(0.5, -0.3)),sd=1.157)
  muestras2=rbind(muestras2,arma_sim)
}

m50=1/500*rowSums(muestras2*muestras2)
print(0000000000)
print(m50)

m51=1/500*rowSums(muestras2[,1:499]*muestras2[,2:500])

print(111111111111111)
print(m51)

m52=1/488*rowSums(muestras2[,1:498]*muestras2[,3:500])

print(2222222222222222)
print(m52)

det1=rbind(c(1,2),c(1,2))


muestrascl=cbind(m50,m51,m52)

nrow(muestrascl)



resultado2=c()

for (i in 1:500 ){
  Y=muestrascl[i,]
  
  
  fobj<-function(par){
    a0=par[1]
    a1=par[2]
    a2=par[3]
    
    M1=matrix(0,2,2)
    M2=matrix(0,2,2)
    
    M1[1,1]=a0^2+a1^2+a2^2
    M1[2,2]=a0^2+a1^2+a2^2
    M1[1,2]=a0*a1+a1*a2
    M1[2,1]=a0*a1+a1*a2
    
    M2[1,1]=a0^2+a1^2+a2^2
    M2[2,2]=a0^2+a1^2+a2^2
    M2[1,2]=a0*a2
    M2[2,1]=a0*a2
    
    S1=matrix(0,2,2)
    S2=matrix(0,2,2)

    
    S1[1,1]=Y[1]
    S1[2,2]=Y[1]
    
    S1[1,2]=Y[2]
    S1[2,1]=Y[2]
    
    S2[1,1]=Y[1]
    S2[2,2]=Y[1]
    
    S2[1,2]=Y[3]
    S2[2,1]=Y[3]
    
    M1inv=chol2inv(chol(M1))
    M2inv=chol2inv(chol(M2))
    
    cl=1/2*log(((a0^2+a1^2+a2^2)^2-(a0*a1+a1*a2)^2))+1/2*sum(diag(M1inv%*%S1))+1/2*log(((a0^2+a1^2+a2^2)^2-(a0*a2)^2))+1/2*sum(diag(M2inv%*%S2))
    return(cl)
  }
  re=optim(c(sd(Y),0,0),fobj)
  
  
  a0=re$par[1]
  a1=re$par[2]
  a2=re$par[3]
  
  resultado2=rbind(resultado2,c(a0,a1,a2))
}
  
meana0=sum(resultado2[,1])/500
meana1=sum(resultado2[,2])/500
meana2=sum(resultado2[,3])/500

biasa0=meana0-1
biasa1=meana1-0.5
biasa2=meana2+0.3

std0=sqrt(sum((resultado2[,1]-meana0)^2)/500)
std1=sqrt(sum((resultado2[,2]-meana1)^2)/500)
std2=sqrt(sum((resultado2[,3]-meana2)^2)/500)

meana0
meana1
meana2



biasa0
biasa1
biasa2

std0  
std1
std2

meana0
meana1
meana2

biasa0
biasa1
biasa2

std0
std1
std2


meana2


#re=optim(c(sd(Y),0,0),fobj)
  
  
#  a0=re$par[1]
#  a1=re$par[2]
#  a2=re$par[3]
  
#  resultado=rbind(resultado,c(a0,a1,a2))
#}









