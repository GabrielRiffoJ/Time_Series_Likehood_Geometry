using HomotopyContinuation
using Symbolics
using LinearAlgebra
using Distributions

@variables a_0 a_1 a_2

n=3

M = zeros(Num, n, n)
    

for i in 1:n
    M[i, i] = a_0
end
    

for i in 2:n
    M[i, i-1] = a_1
    M[i-1, i] = a_1
end
    

for i in 3:n
    M[i, i-2] = a_2
    M[i-2, i] = a_2
end


M0 = zeros(Num, n, n)
    

for i in 1:n
    M0[i, i] = 1
end


M1 = zeros(Num, n, n)


for i in 2:n
    M1[i, i-1] = 1
    M1[i-1, i] = 1
end


M2 = zeros(Num, n, n)


for i in 3:n
    M2[i, i-2] = 1
    M2[i-2, i] = 1
end

using Random, Primes


primes_list = primes(2, 2000)  
primes_list = primes_list[1:300]


reciprocals = 1.0 ./ primes_list


m = 300
e = randn(m)
theta1, theta2 = 0.5, 0.3
ma2 = e .+ theta1 .* [0.0; e[1:end-1]] .+ theta2 .* [0.0; 0.0; e[1:end-2]]


primes_list, reciprocals, ma2


Y=primes_list[1:n]

deter=det(M)

Minv=inv(M)

Minv*M

Minvdet=Minv*deter


l0= tr(Minvdet*M0)*deter-transpose(Y)*Minvdet*M0*Minvdet*Y
l1= tr(Minvdet*M1)*deter-transpose(Y)*Minvdet*M1*Minvdet*Y
l2= tr(Minvdet*M2)*deter-transpose(Y)*Minvdet*M2*Minvdet*Y


f0=simplify(l0,expand=true)
f1=simplify(l1,expand=true)
f2=simplify(l2,expand=true)



@var a0 a1 a2

fa0 = string(f0)
fa1 = string(f1)
fa2 = string(f2)
#fa2 = string(fa1ex)
#fa2_str = string(fa2ex)
#fa3_str = string(fa3ex)



fa0_str = replace(fa0, "a_0" => "a0","a_1" => "a1","a_2" => "a2")
fa1_str = replace(fa1, "a_0" => "a0","a_1" => "a1","a_2" => "a2")
fa2_str = replace(fa2, "a_0" => "a0","a_1" => "a1","a_2" => "a2")



fa0_expr = eval(Meta.parse(fa0_str))
fa1_expr = eval(Meta.parse(fa1_str))
fa2_expr = eval(Meta.parse(fa2_str))


F = System([fa0_expr,fa1_expr,fa2_expr])

t1 = time()

resultcl2 = solve(F)

t2 = time()

elapsed_time = t2 - t1

solutions(resultcl2)

elapsed_time


