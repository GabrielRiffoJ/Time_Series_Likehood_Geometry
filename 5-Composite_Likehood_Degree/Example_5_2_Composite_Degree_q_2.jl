using HomotopyContinuation
using Symbolics
using LinearAlgebra

@variables a_0 a_1 a_2


M1 = zeros(Num, 2, 2)

M1=[ a_0^2+a_1^2+a_2^2   a_0*a_1+a_1*a_2
     a_0*a_1+a_1*a_2   a_0^2+a_1^2+a_2^2
]

M2 = zeros(Num, 2, 2)


M2=[ a_0^2+a_1^2+a_2^2   a_0*a_2
     a_0*a_2   a_0^2+a_1^2+a_2^2
]



@variables x y  
a0M1 = Matrix{Num}(undef, 2, 2)  
fill!(a0M1, 0)  

a1M1 = Matrix{Num}(undef, 2, 2)  
fill!(a1M1, 0)  

a2M1 = Matrix{Num}(undef, 2, 2)  
fill!(a2M1, 0)  




D = a_0^2 + a_1^2 + a_2^2   
U1 = a_0 * a_1 + a_1 * a_2  


a0M1[1,1] = 2*a_0
a0M1[2,2] = 2*a_0

a0M1[1,2] = a_1
a0M1[2,1] = a_1

a1M1[1,1] = 2*a_1
a1M1[2,2] = 2*a_1

a1M1[1,2] = a_0+a_2
a1M1[2,1] = a_0+a_2

a2M1[1,1] = 2*a_2
a2M1[2,2] = 2*a_2

a2M1[1,2] = a_1
a2M1[2,1] = a_1


a0M2 = Matrix{Num}(undef, 2, 2)  
fill!(a0M2, 0) 

a1M2 = Matrix{Num}(undef, 2, 2)  
fill!(a1M2, 0) 

a2M2 = Matrix{Num}(undef, 2, 2)  
fill!(a2M2, 0)  




D = a_0^2 + a_1^2 + a_2^2   
U2 = a_0 * a_2   


a0M2[1,1] = 2*a_0
a0M2[2,2] = 2*a_0

a0M2[1,2] = a_2
a0M2[2,1] = a_2

a1M2[1,1] = 2*a_1
a1M2[2,2] = 2*a_1

a1M2[1,2] = 0
a1M2[2,1] = 0

a2M2[1,1] = 2*a_2
a2M2[2,2] = 2*a_2

a2M2[1,2] = a_0
a2M2[2,1] = a_0

SM0=[1   0
     0   1
]

SMH=[0   1
     1   0
]

S1=[ 
1/3 1/13
1/13 1/3
]

S2=[
1/3 1/17
1/17 1/3
]

M1
M2

deterM=1
deterM=deterM*det(M1)*det(M2)
M1inv=inv(M1)*det(M1)
M2inv=inv(M2)*det(M2)

detM1=det(M1)
detM1ex=Symbolics.simplify(detM1,expand=true)

detM2=det(M2)
detM2ex=Symbolics.simplify(detM2,expand=true)



f1a0 = detM1*(-0.5*tr(a0M1*(M1inv))) -0.5* tr( (-(M1inv) * a0M1 * (M1inv)) * S1) 
f1a1 = detM1*(-0.5*tr(a1M1*(M1inv))) -0.5* tr( (-(M1inv) * a1M1 * (M1inv)) * S1)
f1a2 = detM1*(-0.5*tr(a2M1*(M1inv))) -0.5* tr( (-(M1inv) * a2M1 * (M1inv)) * S1)

f1a0ex=Symbolics.simplify(detM2^2*f1a0,expand=true)
f1a1ex=Symbolics.simplify(detM2^2*f1a1,expand=true)
f1a2ex=Symbolics.simplify(detM2^2*f1a2,expand=true)


f2a0 = (detM2)*(n-2)*(-0.5*tr(a0M2*(M2inv))) -0.5* tr( (-(M2inv) * a0M2 * (M2inv)) * S2) 
f2a1 = (detM2)*(n-2)*(-0.5*tr(a1M2*(M2inv))) -0.5* tr( (-(M2inv) * a1M2 * (M2inv)) * S2)
f2a2 = (detM2)*(n-2)*(-0.5*tr(a2M2*(M2inv))) -0.5* tr( (-(M2inv) * a2M2 * (M2inv)) * S2)

f2a0ex=Symbolics.simplify(detM1^2*f2a0,expand=true)
f2a1ex=Symbolics.simplify(detM1^2*f2a1,expand=true)
f2a2ex=Symbolics.simplify(detM1^2*f2a2,expand=true)


fa0=f1a0ex+f2a0ex
fa1=f1a1ex+f2a1ex
fa2=f1a2ex+f2a2ex


fa0ex=Symbolics.simplify(fa0,expand=true)
fa1ex=Symbolics.simplify(fa1,expand=true)
fa2ex=Symbolics.simplify(fa2,expand=true)




@var a0 a1 a2 

fa0 = string(fa0ex)
fa1 = string(fa1ex)
fa2 = string(fa2ex)


fa0_str = replace(fa0, "a_0" => "a0","a_1" => "a1","a_2" => "a2")
fa1_str = replace(fa1, "a_0" => "a0","a_1" => "a1","a_2" => "a2")
fa2_str = replace(fa2, "a_0" => "a0","a_1" => "a1","a_2" => "a2")

fa0_expr = eval(Meta.parse(fa0_str))
fa1_expr = eval(Meta.parse(fa1_str))
fa2_expr = eval(Meta.parse(fa2_str))

F = System([fa0_expr,fa1_expr,fa2_expr])

t1 = Sys.time()

resultcla_MA2 = HomotopyContinuation.solve(F)

t2=Sys.time()


#(t2-t1)/3600*500
#nonsingular(resultcla_MA2)


#resultcla_MA2


nonsingular(resultcla_MA2)

noninvertible1=[]
noninvertible1_gamma=[]
noninvertible2=[]
noninvertible2_gamma=[]
noninvertible3=[]
noninvertible3_gamma=[]
noninvertiblea_0=[]
noninvertiblea_0_gamma=[]
invertible=[]
invertible_gamma=[]

totsol=length(nonsingular(resultcla_MA2))
tol=10^(-6)

i=1
for i in 1:totsol
     solutiontemp=nonsingular(resultcla_MA2)[i].solution
     a_0=solutiontemp[1]
     a_1=solutiontemp[2]
     a_2=solutiontemp[3]
     gamma0=a_0^2+a_1^2+a_2^2
     gamma1=a_0*a_1+a_1*a_2
     gamma2=a_0*a_2

     if abs(a_1)<tol
        if abs(a_0+a_2)<tol
            push!(noninvertiblea_0,[a_0,a_1,a_2])
            push!(noninvertiblea_0_gamma,[gamma0,gamma1,gamma2])
        end
    elseif abs(a_0+a_1+a_2)<tol
        push!(noninvertible1,[a_0,a_1,a_2])
        push!(noninvertible1_gamma,[gamma0,gamma1,gamma2])
    elseif abs(a_0-a_1+a_2)<tol
        push!(noninvertible2,[a_0,a_1,a_2])
        push!(noninvertible2_gamma,[gamma0,gamma1,gamma2])
    elseif abs(a_0-a_2)<tol
        push!(noninvertible3,[a_0,a_1,a_2])
        push!(noninvertible3_gamma,[gamma0,gamma1,gamma2])
    else 
        push!(invertible,[a_0,a_1,a_2])
        push!(invertible_gamma,[gamma0,gamma1,gamma2])
    end
end

noninvertible1
noninvertible2
noninvertible3
noninvertiblea_0
invertible


using LinearAlgebra

function unique_vectors_tol(vectors; tol = 1e-8)
    uniques = Vector{eltype(vectors)}()
    for v in vectors
         if all(norm(v - u) > tol for u in uniques)
            push!(uniques, v)
        end
    end
    return uniques
end


noninvertible1_gamma_clean   = unique_vectors_tol(noninvertible1_gamma; tol=tol)
noninvertible2_gamma_clean   = unique_vectors_tol(noninvertible2_gamma; tol=tol)
noninvertible3_gamma_clean   = unique_vectors_tol(noninvertible3_gamma; tol=tol)
noninvertiblea_0_gamma_clean = unique_vectors_tol(noninvertiblea_0_gamma; tol=tol)
invertible_gamma_clean       = unique_vectors_tol(invertible_gamma; tol=tol)

noninvertible1_gamma_clean  
noninvertible2_gamma_clean   
noninvertible3_gamma_clean   
noninvertiblea_0_gamma_clean 
invertible_gamma_clean       



gamma_1