using HomotopyContinuation
using Symbolics
using LinearAlgebra



@variables a_1 a_2 a_3 a_4 a_5 a_6 a_7 a_8 a_9 a_10

n=20

M = zeros(Num, n, n)
    

for i in 1:n
    M[i, i] = a_2^2+a_1^2+1
end

for i in 1:(n-1)
    M[i, i+1] = -a_1+a_1*a_2
    M[i+1, i] = -a_1+a_1*a_2
end

for i in 1:(n-2)
    M[i, i+2] = -a_2
    M[i+2, i] = -a_2
end


D1 = zeros(Num, n, n)

for i in 1:(n)
    D1[i, i] = 2*a_1
end


for i in 1:(n-1)
    D1[i, i+1] = -1+a_2
    D1[i+1, i] = -1+a_2
end

D2 = zeros(Num, n, n)

for i in 1:(n)
    D2[i, i] = 2*a_2
end


for i in 1:(n-1)
    D2[i, i+1] = a_1
    D2[i+1, i] = a_1
end

for i in 1:(n-2)
    D2[i, i+2] = -1
    D2[i+2, i] = -1
end

using Random
using Primes


primes = collect(primes(1, 2000))      
primes_300 = primes[1:300]             

inv_primes = 1.0 ./ primes_300


m = 300               
phi1 = 0.5            
phi2 = 0.3            
sigma = 1.0           

ar2 = zeros(m)
epsilon = sigma .* randn(n)


ar2[1] = epsilon[1]
ar2[2] = phi1*ar2[1] + epsilon[2]


for t in 3:m
    ar2[t] = phi1*ar2[t-1] + phi2*ar2[t-2] + epsilon[t]
end

#Y=primes_300[1:n]

Y=rand(-100:100,n)
#Y=primes_300[rand(1:300,n)]
#Y=primes_300[rand(1:300,n)]



M[1,1]=1
M[1,2]=-a_1
M[2,1]=-a_1
M[2,2]=a_1^2+1

M[n,n]=1
M[n-1,n]=-a_1
M[n,n-1]=-a_1
M[n-1,n-1]=a_1^2+1

D1[1,1]=0
D1[1,2]=-1
D1[2,1]=-1
D1[2,2]=2*a_1

D1[n,n]=0
D1[n-1,n]=-1
D1[n,n-1]=-1
D1[n-1,n-1]=2*a_1


D2[1,1]=0
D2[1,2]=0
D2[2,1]=0
D2[2,2]=0

D2[n,n]=0
D2[n-1,n]=0
D2[n,n-1]=0
D2[n-1,n-1]=0

Minv=inv(M)
deter=det(M)
Minvdet=Minv*deter

for i in 1:n
    for j in 1:n
        Minvdet[i,j]=simplify(Minv[i,j]*deter,expand=true)
    end
end


@variables x y
@variables a_1 a_2
    
l1= tr(D1*Minvdet)*transpose(Y)*M*Y-n*transpose(Y)*D1*Y*deter
l2= tr(D2*Minvdet)*transpose(Y)*M*Y-n*transpose(Y)*D2*Y*deter

f1=simplify(l1,expand=true)
f2=simplify(l2,expand=true)

fa1 = string(f1)
fa2 = string(f2)

@var a0 a1 a2

fa1_str = replace(fa1,"a_1" => "a1","a_2" => "a2")
fa2_str = replace(fa2,"a_1" => "a1","a_2" => "a2")

fa1_expr = eval(Meta.parse(fa1_str))
fa2_expr = eval(Meta.parse(fa2_str))

F = System([fa1_expr,fa2_expr])

resultcl2 = HomotopyContinuation.solve(F)


