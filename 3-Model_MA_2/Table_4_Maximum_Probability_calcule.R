n <- 75

n1=ceiling(n/2)
n2=floor(n/2) 

R=matrix(0,n,n)
diag(R)=1
for (i in 1:(n-2)){
  R[i,i+2]=-1/2
  R[i+2,i]=-1/2
}

Rinv=solve(R)


R11=matrix(0,n,n)
diag(R11)=0
for (i in 1:(n-2)){
  R11[i,i+2]=1/2
  R11[i+2,i]=1/2
}


R12=matrix(0,n,n)
diag(R12)=0
for (i in 1:(n-1)){
  R12[i,i+1]=1/2
  R12[i+1,i]=1/2
}

R2=matrix(0,n,n)
diag(R2)=0
for (i in 1:(n-2)){
    R2[i,i+2]=-1/2
    R2[i+2,i]=-1/2
}

R22=matrix(0,n,n)
diag(R22)=0
for (i in 1:(n-2)){
  R22[i,i+2]=1/2
  R22[i+2,i]=1/2
}


#H11
H11t=-sum(Rinv*R11)
H11m=Rinv%*%R11%*%Rinv

H12t=-sum(Rinv*R12)
H12m= Rinv%*%R12%*%Rinv



H22t=-sum(diag((-Rinv%*%R2%*%Rinv%*%R2+Rinv%*%R22)))
H22m1=Rinv%*%R2%*%Rinv%*%R2%*%Rinv
H22m2=Rinv%*%R22%*%Rinv
H22m3=Rinv%*%R2%*%Rinv


H75=c()



theta1 <- 0
theta2 <- -1

v1=c()
v2=c()



for (i in 1:10000){



  Y <- arima.sim(model = list(ma = c(theta1, theta2)), n = n)

  H11=H11t+n*(t(Y)%*%H11m%*%Y)/(t(Y)%*%Rinv%*%Y)
  H12=H12t+n*(t(Y)%*%H12m%*%Y)/(t(Y)%*%Rinv%*%Y)
  H22=H12t+n*1/(t(Y)%*%Rinv%*%Y)*(-2*t(Y)%*%H22m1%*%Y*t(Y)%*%Rinv%*%Y+t(Y)%*%H22m2%*%Y*t(Y)%*%Rinv%*%Y+(t(Y)%*%H22m3%*%Y)^2)

  v1=c(v1,H11)
  v2=c(v2,H11*H22-(H12)^2)

}

hist(v1)
hist(v2)
  
sum(v1<0)
sum(v2>0)

H75=cbind(H75,v1,v2)


colSums(H<0)
colSums(H25<0)
colSums(H50<0)
colSums(H75<0)
colSums(H100<0)
colSums(H500<0)
H25[1,]


