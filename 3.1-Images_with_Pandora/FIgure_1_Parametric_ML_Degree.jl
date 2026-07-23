using HomotopyContinuation
using Symbolics

@variables a_0 a_1 a_2

using Symbolics




n=5
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
U2a2 = a_2  
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


Y=[2,3,5,7,11,13,17,19,23,29,31,37]

@variables y_[1:10]



Y_vec=Y[1:n]
Y_vec=[y_[1],y_[2],y_[3],y_[4],y_[5]]



M
Ma0
Ma1
Ma2

detM=det(M)
detMex=Symbolics.simplify(detM,expand=true)

Minv2=(detM*Minv)
Minv2_sim=Minv2
for i in 1:n
    for j in 1:n
        Minv2_sim[i,j] = Symbolics.simplify(Minv2[i,j],expand=true)
    end
end
Minv2ex=Minv2_sim



tra0=detMex*(-0.5*tr(Ma0*Minv2ex))
tra1=detMex*(-0.5*tr(Ma1*Minv2ex))
tra2=detMex*(-0.5*tr(Ma2*Minv2ex))

traa0_sim=Symbolics.simplify(tra0,expand=true)
traa1_sim=Symbolics.simplify(tra1,expand=true)
traa2_sim=Symbolics.simplify(tra2,expand=true)

bila0= -0.5* transpose(Y_vec) * (-Minv2ex * Ma0*Minv2ex) * Y_vec  
bila1= -0.5* transpose(Y_vec) * (-Minv2ex * Ma1*Minv2ex) * Y_vec
bila2= -0.5* transpose(Y_vec) * (-Minv2ex * Ma2*Minv2ex) * Y_vec

bila0_sim= Symbolics.simplify(bila0,expand=true)
bila1_sim= Symbolics.simplify(bila1,expand=true)
bila2_sim= Symbolics.simplify(bila2,expand=true)

fa0_temp=traa0_sim+bila0_sim
fa1_temp=traa1_sim+bila1_sim
fa2_temp=traa2_sim+bila2_sim

fa0_sim=Symbolics.simplify(fa0_temp,expand=true)
fa1_sim=Symbolics.simplify(fa1_temp,expand=true)
fa2_sim=Symbolics.simplify(fa2_temp,expand=true)


fa0 = detMex*(-0.5*tr(Ma0*Minv2ex)) -0.5* transpose(Y_vec) * (-Minv2ex * Ma0*Minv2ex) * Y_vec  
fa1 = detMex*(-0.5*tr(Ma1*Minv2ex)) -0.5* transpose(Y_vec) * (-Minv2ex * Ma1*Minv2ex) * Y_vec
fa2 = detMex*(-0.5*tr(Ma2*Minv2ex)) -0.5* transpose(Y_vec) * (-Minv2ex * Ma2*Minv2ex) * Y_vec



fa0ex2=Symbolics.simplify(fa0,expand=true)
fa1ex2=Symbolics.simplify(fa1,expand=true)
fa2ex2=Symbolics.simplify(fa2,expand=true)

#fa0ex2=fa0_sim
#fa1ex2=fa1_sim
#fa2ex2=fa2_sim

using HomotopyContinuation
@var a0 a1 a2

fa0_str = string(fa0ex2)
fa1_str = string(fa1ex2)
fa2_str = string(fa2ex2)


@var a0 a1 a2 y[1:10]


fa0_str = replace(fa0_str, "a_0" => "a0", "a_1" => "a1", "a_2" => "a2", "y_[1]" => "y[1]", "y_[2]" => "y[2]", "y_[3]" => "y[3]", "y_[4]" => "y[4]", "y_[5]" => "y[5]", "y_[6]" => "y[6]", "y_[7]" => "y[7]", "y_[8]" => "y[8]", "y_[9]" => "y[9]", "y_[10]" => "y[10]")
fa1_str = replace(fa1_str, "a_0" => "a0", "a_1" => "a1", "a_2" => "a2", "y_[1]" => "y[1]", "y_[2]" => "y[2]", "y_[3]" => "y[3]", "y_[4]" => "y[4]", "y_[5]" => "y[5]", "y_[6]" => "y[6]", "y_[7]" => "y[7]", "y_[8]" => "y[8]", "y_[9]" => "y[9]", "y_[10]" => "y[10]")
fa2_str = replace(fa2_str, "a_0" => "a0", "a_1" => "a1", "a_2" => "a2", "y_[1]" => "y[1]", "y_[2]" => "y[2]", "y_[3]" => "y[3]", "y_[4]" => "y[4]", "y_[5]" => "y[5]", "y_[6]" => "y[6]", "y_[7]" => "y[7]", "y_[8]" => "y[8]", "y_[9]" => "y[9]", "y_[10]" => "y[10]")



fa0_expr = eval(Meta.parse(fa0_str))
fa1_expr = eval(Meta.parse(fa1_str))
fa2_expr = eval(Meta.parse(fa2_str))



F = System([fa0_expr, fa1_expr,fa2_expr])



t1 = time()


result12 = HomotopyContinuation.solve(F)

t2 = time()
elapsed_time = t2 - t1

YTime[n]=elapsed_time



using Pandora

E = EnumerativeProblem([fa0_expr,fa1_expr, fa2_expr], variables=[a0,a1, a2], parameters=vcat(y), torus_only=true)

ER = tally.(explore(E, [n_real_solutions, n_positive_solutions], n_samples=1000))
O = maximize_n_real_solutions(E)
C = certify(record_fibre(O), E)
(V, P) = Pandora.visualize(E; near=record_parameters(O), strategy=:quadtree)

refine!(V)
refine!(V)
P2 = Pandora.visualize(V)



using Plots: savefig



using Plots


#savefig(P2, "MA2_an3.png")
#Plots.savefig(P2, "MA2_an3.png")

display(P2)



