using HomotopyContinuation
using Symbolics
@variables a_0 a_1 a_2 a_3
using Pkg
using Distributions

Pkg.add("Primes")
using Primes

prime = collect(primes(1, 1987))[1:300]

#YTime = Vector{Float64}(undef, 100)

theta1 = 0.6
theta2 = -0.2
theta3 = -0.1
sigma  = 1.0

m = 300


ep = rand(Normal(0, sigma), m + 2)


Ym = zeros(m)
for t in 1:m
    Ym[t] = ep[t+2] + theta1 * ep[t+1] + theta2 * ep[t]
end


n=4

YI=prime
YF=  1.0 ./ prime

#first approach
Y=YF[1:n]

@variables x y 
M = Matrix{Num}(undef, n, n)  
fill!(M, 0)  


D  = a_0^2 + a_1^2 + a_2^2 + a_3^2
U1 = a_0*a_1 + a_1*a_2 + a_2*a_3
U2 = a_0*a_2 + a_1*a_3
U3 = a_0*a_3

for i in 1:n
    M[i, i] = D

    if i < n
        M[i, i+1] = U1
        M[i+1, i] = U1
    end

    if i < n-1
        M[i, i+2] = U2
        M[i+2, i] = U2
    end

    if i < n-2
        M[i, i+3] = U3
        M[i+3, i] = U3
    end
end


Minv=inv(M)

Ma0 = Matrix{Num}(undef, n, n)
fill!(Ma0, 0)

Da0  = 2*a_0
U1a0 = a_1
U2a0 = a_2
U3a0 = a_3

for i in 1:n
    Ma0[i,i] = Da0

    if i < n
        Ma0[i,i+1] = U1a0
        Ma0[i+1,i] = U1a0
    end

    if i < n-1
        Ma0[i,i+2] = U2a0
        Ma0[i+2,i] = U2a0
    end

    if i < n-2
        Ma0[i,i+3] = U3a0
        Ma0[i+3,i] = U3a0
    end
end


Ma1 = Matrix{Num}(undef, n, n)
fill!(Ma1, 0)

Da1  = 2*a_1
U1a1 = a_0 + a_2
U2a1 = a_3
U3a1 = 0

for i in 1:n
    Ma1[i,i] = Da1

    if i < n
        Ma1[i,i+1] = U1a1
        Ma1[i+1,i] = U1a1
    end

    if i < n-1
        Ma1[i,i+2] = U2a1
        Ma1[i+2,i] = U2a1
    end

    if i < n-2
        Ma1[i,i+3] = U3a1
        Ma1[i+3,i] = U3a1
    end
end


Ma2 = Matrix{Num}(undef, n, n)
fill!(Ma2, 0)

Da2  = 2*a_2
U1a2 = a_1 + a_3
U2a2 = a_0
U3a2 = 0

for i in 1:n
    Ma2[i,i] = Da2

    if i < n
        Ma2[i,i+1] = U1a2
        Ma2[i+1,i] = U1a2
    end

    if i < n-1
        Ma2[i,i+2] = U2a2
        Ma2[i+2,i] = U2a2
    end

    if i < n-2
        Ma2[i,i+3] = U3a2
        Ma2[i+3,i] = U3a2
    end
end


Ma3 = Matrix{Num}(undef, n, n)
fill!(Ma3, 0)

Da3  = 2*a_3
U1a3 = a_2
U2a3 = a_1
U3a3 = a_0

for i in 1:n
    Ma3[i,i] = Da3

    if i < n
        Ma3[i,i+1] = U1a3
        Ma3[i+1,i] = U1a3
    end

    if i < n-1
        Ma3[i,i+2] = U2a3
        Ma3[i+2,i] = U2a3
    end

    if i < n-2
        Ma3[i,i+3] = U3a3
        Ma3[i+3,i] = U3a3
    end
end


using LinearAlgebra
using Symbolics

#import Pkg
#Pkg.add("SymPy")

#using SymPy

#@variables Y[1:n] 

Y_vec=Y[1:n]


detM=det(M)
detMex=simplify(detM,expand=true)

Minv3=(detM*Minv)
Minv3ex=Minv3




fa0 = detMex*(-0.5*tr(Ma0*Minv3ex)) -0.5* (transpose(Y)*(-Minv3ex * Ma0*Minv3ex) * Y)  
fa1 = detMex*(-0.5*tr(Ma1*Minv3ex)) -0.5* (transpose(Y)*(-Minv3ex * Ma1*Minv3ex) * Y)
fa2 = detMex*(-0.5*tr(Ma2*Minv3ex)) -0.5* (transpose(Y)*(-Minv3ex * Ma2*Minv3ex) * Y)
fa3 = detMex*(-0.5*tr(Ma3*Minv3ex)) -0.5* (transpose(Y)*(-Minv3ex * Ma3*Minv3ex) * Y)

#fa0 = detMex*(-0.5*tr(Ma0*Minv2ex)) -0.5* tr(transpose(Y)*(-Minv2ex * Ma0*Minv2ex) * Y)  
#fa1 = detMex*(-0.5*tr(Ma1*Minv2ex)) -0.5* tr(transpose(Y)*(-Minv2ex * Ma1*Minv2ex) * Y)
#fa2 = detMex*(-0.5*tr(Ma2*Minv2ex)) -0.5* tr(transpose(Y)*(-Minv2ex * Ma2*Minv2ex) * Y)

fa0ex2=simplify(fa0,expand=true)
fa1ex2=simplify(fa1,expand=true)
fa2ex2=simplify(fa2,expand=true)
fa3ex2=simplify(fa3,expand=true)

using HomotopyContinuation
@var a0 a1 a2 a3

fa0_str = string(fa0ex2)
fa1_str = string(fa1ex2)
fa2_str = string(fa2ex2)
fa3_str = string(fa3ex2)


@var a0 a1 a2  a3


fa0_str = replace(fa0_str, "a_0" => "a0", "a_1" => "a1", "a_2" => "a2", "a_3" => "a3")
fa1_str = replace(fa1_str, "a_0" => "a0", "a_1" => "a1", "a_2" => "a2", "a_3" => "a3")
fa2_str = replace(fa2_str, "a_0" => "a0", "a_1" => "a1", "a_2" => "a2", "a_3" => "a3")
fa3_str = replace(fa3_str, "a_0" => "a0", "a_1" => "a1", "a_2" => "a2", "a_3" => "a3")


fa0_expr = eval(Meta.parse(fa0_str))
fa1_expr = eval(Meta.parse(fa1_str))
fa2_expr = eval(Meta.parse(fa2_str))
fa3_expr = eval(Meta.parse(fa3_str))


F = System([fa0_expr, fa1_expr,fa2_expr, fa3_expr])



t1 = time()


resultS = solve(F)

t2 = time()

elapsed_time = t2 - t1

elapsed_time

