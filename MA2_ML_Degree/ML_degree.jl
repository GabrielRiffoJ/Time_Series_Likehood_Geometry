using HomotopyContinuation
using Symbolics
@variables a_0 a_1 a_2
using Pkg

Pkg.add("Primes")
using Primes

prime = collect(primes(1, 1987))[1:300]

#YTime = Vector{Float64}(undef, 100)

theta1 = 0.6
theta2 = -0.3
sigma  = 1.0

m = 300


ep = rand(Normal(0, sigma), m + 2)


Ym = zeros(m)
for t in 1:m
    Ym[t] = ep[t+2] + theta1 * ep[t+1] + theta2 * ep[t]
end


n=3

YI=prime
YF=  1.0 ./ prime

#first approach
Y=YI[1:n]

@variables x y 
M = Matrix{Num}(undef, n, n)  
fill!(M, 0)  


D = a_0^2 + a_1^2 + a_2^2  
U1 = a_0 * a_1 + a_1 * a_2   
U2 = a_0 * a_2  

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
end



Minv=inv(M)

Ma0 = Matrix{Num}(undef, n, n)  
fill!(Ma0, 0)  
Da0 = 2a_0 
U1a0 = a_1 
U2a0 =  a_2 
for i in 1:n
    Ma0[i, i] = Da0  
    if i < n
        Ma0[i, i+1] = U1a0  
        Ma0[i+1, i] = U1a0  
    end
    if i < n-1
        Ma0[i, i+2] = U2a0  
        Ma0[i+2, i] = U2a0  
    end
end



Ma1 = Matrix{Num}(undef, n, n) 
fill!(Ma1, 0) 
Da1 = 2*a_1   
U1a1 = a_0+ a_2 
U2a1 = 0  
for i in 1:n
    Ma1[i, i] = Da1 
    if i < n
        Ma1[i, i+1] = U1a1 
        Ma1[i+1, i] = U1a1 
    end
    if i < n-1
        Ma1[i, i+2] = U2a1 
        Ma1[i+2, i] = U2a1 
    end
end



Ma2 = Matrix{Num}(undef, n, n) 
fill!(Ma2, 0) 
Da2 =  2*a_2  
U1a2 =  a_1  
U2a2 = a_0 
for i in 1:n
    Ma2[i, i] = Da2  
    if i < n
        Ma2[i, i+1] = U1a2  
        Ma2[i+1, i] = U1a2  
    end
    if i < n-1
        Ma2[i, i+2] = U2a2  
        Ma2[i+2, i] = U2a2  
    end
end


using LinearAlgebra
using Symbolics

#import Pkg
#Pkg.add("SymPy")

#using SymPy

#@variables Y[1:n]  # Definir un vector de n parámetros

Y_vec=Y[1:n]




detM=det(M)
detMex=simplify(detM,expand=true)

Minv2=(detM*Minv)
Minv2ex=Minv2




fa0 = detMex*(-0.5*tr(Ma0*Minv2ex)) -0.5* (transpose(Y)*(-Minv2ex * Ma0*Minv2ex) * Y)  
fa1 = detMex*(-0.5*tr(Ma1*Minv2ex)) -0.5* (transpose(Y)*(-Minv2ex * Ma1*Minv2ex) * Y)
fa2 = detMex*(-0.5*tr(Ma2*Minv2ex)) -0.5* (transpose(Y)*(-Minv2ex * Ma2*Minv2ex) * Y)


#fa0 = detMex*(-0.5*tr(Ma0*Minv2ex)) -0.5* tr(transpose(Y)*(-Minv2ex * Ma0*Minv2ex) * Y)  
#fa1 = detMex*(-0.5*tr(Ma1*Minv2ex)) -0.5* tr(transpose(Y)*(-Minv2ex * Ma1*Minv2ex) * Y)
#fa2 = detMex*(-0.5*tr(Ma2*Minv2ex)) -0.5* tr(transpose(Y)*(-Minv2ex * Ma2*Minv2ex) * Y)

fa0ex2=simplify(fa0,expand=true)
fa1ex2=simplify(fa1,expand=true)
fa2ex2=simplify(fa2,expand=true)


using HomotopyContinuation
@var a0 a1 a2

fa0_str = string(fa0ex2)
fa1_str = string(fa1ex2)
fa2_str = string(fa2ex2)


@var a0 a1 a2  


fa0_str = replace(fa0_str, "a_0" => "a0", "a_1" => "a1", "a_2" => "a2")
fa1_str = replace(fa1_str, "a_0" => "a0", "a_1" => "a1", "a_2" => "a2")
fa2_str = replace(fa2_str, "a_0" => "a0", "a_1" => "a1", "a_2" => "a2")



fa0_expr = eval(Meta.parse(fa0_str))
fa1_expr = eval(Meta.parse(fa1_str))
fa2_expr = eval(Meta.parse(fa2_str))



F = System([fa0_expr, fa1_expr,fa2_expr])



t1 = time()


resultS = solve(F)

t2 = time()

elapsed_time = t2 - t1

elapsed_time

for i in 1:58
    println(solutions(resultS)[i])
end

for i in 1:58
    solg1=solutions(resultS)[i][2]
    if abs(solg1)<0.001
        println(i)
    end
end

contsim1=0
for i in 1:58
    solg0=solutions(resultS)[i][1]
    solg1=solutions(resultS)[i][2]
    solg2=solutions(resultS)[i][3]
    if abs(solg0-solg2)<0.001
        contsim1=contsim1+1
        println(i)
    end
end

contsim1


contsim4=0
for i in 1:58
    solg0=solutions(resultS)[i][1]
    solg1=solutions(resultS)[i][2]
    solg2=solutions(resultS)[i][3]
    if abs(solg0+solg2)<0.001
        contsim4=contsim4+1
        println(i)
    end
end

contsim4

contsim2=0
for i in 1:58
    solg0=solutions(resultS)[i][1]
    solg1=solutions(resultS)[i][2]
    solg2=solutions(resultS)[i][3]
    if abs(solg0-solg1+solg2)<0.001
        contsim2=contsim2+1
        println(i)
    end
end
contsim2

contsim3=0
for i in 1:58
    solg0=solutions(resultS)[i][1]
    solg1=solutions(resultS)[i][2]
    solg2=solutions(resultS)[i][3]
    if abs(solg0+solg1+solg2)<0.001
        contsim3=contsim3+1
        println(i)
    end
end

println(contsim1)
println(contsim2)
println(contsim3)
println(contsim4)


