using HomotopyContinuation
using Symbolics
using LinearAlgebra


@variables a_1 a_2 a_3 a_4 a_5 a_6 a_7 a_8 a_9 a_10

n=21

M = zeros(Num, n, n)
    

for i in 1:n
    M[i, i] = a_3^2+a_2^2+a_1^2+1
end

for i in 1:(n-1)
    M[i, i+1] = -a_1+a_1*a_2+a_2*a_3
    M[i+1, i] = -a_1+a_1*a_2+a_2*a_3
end

for i in 1:(n-2)
    M[i, i+2] = -a_2+a_1*a_3
    M[i+2, i] = -a_2+a_1*a_3
end

for i in 1:(n-3)
    M[i, i+3] = -a_3
    M[i+3, i] = -a_3
end

D1 = zeros(Num, n, n)

for i in 1:(n)
    D1[i, i] = 2*a_1
end


for i in 1:(n-1)
    D1[i, i+1] = -1+a_2
    D1[i+1, i] = -1+a_2
end

for i in 1:(n-2)
    D1[i, i+2] = a_3
    D1[i+2, i] = a_3
end

D2 = zeros(Num, n, n)

for i in 1:(n)
    D2[i, i] = 2*a_2
end


for i in 1:(n-1)
    D2[i, i+1] = a_1+a_3
    D2[i+1, i] = a_1+a_3
end

for i in 1:(n-2)
    D2[i, i+2] = -1
    D2[i+2, i] = -1
end

D3 = zeros(Num, n, n)

for i in 1:(n)
    D3[i, i] = 2*a_3
end


for i in 1:(n-1)
    D3[i, i+1] = a_2
    D3[i+1, i] = a_2
end

for i in 1:(n-2)
    D3[i, i+2] = a_1
    D3[i+2, i] = a_1
end

for i in 1:(n-3)
    D3[i, i+3] = -1
    D3[i+3, i] = -1
end



using Primes
using Random


n_primes = 300
prime_numbers = primes(n_primes)

inv_primes = 1.0 ./ prime_numbers

phi1 = 0.5
phi2 = -0.3
phi3 = 0.2

T = 300                      
n_realizations = 300         


ar3_data = zeros(n_realizations, T)

for r in 1:n_realizations
    ε = randn(T)             
    x = zeros(T)

    for t in 4:T
        x[t] = phi1*x[t-1] + phi2*x[t-2] + phi3*x[t-3] + ε[t]
    end

    ar3_data[r, :] = x
end


Y=YI[1:n]


M[1,1]=1
M[1,2]=-a_1
M[2,1]=-a_1
M[2,2]=a_1^2+1
M[1,3]=-a_2
M[3,1]=-a_2
M[2,3]=-a_1+a_1*a_2
M[3,2]=-a_1+a_1*a_2
M[3,3]=a_1^2+a_2^2+1

M[n,n]=1
M[n,n-1]=-a_1
M[n-1,n]=-a_1
M[n-1,n-1]=a_1^2+1
M[n,n-2]=-a_2
M[n-2,n]=-a_2
M[n-1,n-2]=-a_1+a_1*a_2
M[n-2,n-1]=-a_1+a_1*a_2
M[n-2,n-2]=a_1^2+a_2^2+1


D1[1,1]=0
D1[1,2]=-1
D1[2,1]=-1
D1[2,2]=2*a_1
D1[1,3]=0
D1[3,1]=0
D1[2,3]=-1+a_2
D1[3,2]=-1+a_2
D1[3,3]=2*a_1

D1[n,n]=0
D1[n,n-1]=-1
D1[n-1,n]=-1
D1[n-1,n-1]=2*a_1
D1[n,n-2]=0
D1[n-2,n]=0
D1[n-1,n-2]=-1+a_2
D1[n-2,n-1]=-1+a_2
D1[n-2,n-2]=2*a_1

D2[1,1]=0
D2[1,2]=0
D2[2,1]=0
D2[2,2]=
D2[1,3]=-1
D2[3,1]=-1
D2[2,3]=a_1
D2[3,2]=a_1
D2[3,3]=2*a_2

D2[n,n]=0
D2[n,n-1]=0
D2[n-1,n]=0
D2[n-1,n-1]=0
D2[n,n-2]=-1
D2[n-2,n]=-1
D2[n-1,n-2]=a_1
D2[n-2,n-1]=a_1
D2[n-2,n-2]=2*a_2

D3[1,1]=0
D3[1,2]=0
D3[2,1]=0
D3[2,2]=0
D3[1,3]=0
D3[3,1]=0
D3[2,3]=0
D3[3,2]=0
D3[3,3]=0

D3[n,n]=0
D3[n,n-1]=0
D3[n-1,n]=0
D3[n-1,n-1]=0
D3[n,n-2]=0
D3[n-2,n]=0
D3[n-1,n-2]=0
D3[n-2,n-1]=0
D3[n-2,n-2]=0


deter=gd
Minvdet=Minv

#for i in 1:n
#    for j in 1:n
#        Minvdet[i,j]=simplify(Minv[i,j]*deter,expand=true)
#    end
#end

@variables x y
@variables a_1 a_2
    
l1= tr(D1*Minvdet)*transpose(Y)*M*Y-n*transpose(Y)*D1*Y*deter
l2= tr(D2*Minvdet)*transpose(Y)*M*Y-n*transpose(Y)*D2*Y*deter
l3= tr(D3*Minvdet)*transpose(Y)*M*Y-n*transpose(Y)*D3*Y*deter

f1=simplify(l1,expand=true)
f2=simplify(l2,expand=true)
f3=simplify(l3,expand=true)

fa1 = string(f1)
fa2 = string(f2)
fa3 = string(f3)

@var a0 a1 a2 a3

fa1_str = replace(fa1,"a_1" => "a1","a_2" => "a2","a_3" => "a3")
fa2_str = replace(fa2,"a_1" => "a1","a_2" => "a2","a_3" => "a3")
fa3_str = replace(fa3,"a_1" => "a1","a_2" => "a2","a_3" => "a3")

fa1_expr = eval(Meta.parse(fa1_str))
fa2_expr = eval(Meta.parse(fa2_str))
fa3_expr = eval(Meta.parse(fa3_str))

F = System([fa1_expr,fa2_expr,fa3_expr])

resultar3 = solve(F)

