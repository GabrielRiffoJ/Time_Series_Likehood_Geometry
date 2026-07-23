using HomotopyContinuation
using Symbolics
using LinearAlgebra


@var g0 g1 g2 g3 g4 g5 g6

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


deterM=1
deterM=deterM*det(M1)*det(M2)
M1inv=inv(M1)*det(M1)
M2inv=inv(M2)*det(M2)

detM1=det(M1)
detM1ex=Symbolics.simplify(detM1,expand=true)

detM2=det(M2)
detM2ex=Symbolics.simplify(detM2,expand=true)

deterM=1
deterM=deterM*det(MA1)*det(MA2)
M1inv=inv(MA1)*det(MA1)
M2inv=inv(MA2)*det(MA2)



estim=[]



for act in 1:500
   
    gcl0=muestrasg0[act]
    gcl1=muestrasg1[act]
    gcl2=muestrasg2[act]

    S1=[
        gcl0   gcl1
        gcl1   gcl0
    ]
    S2=[
        gcl0  gcl2
        gcl2  gcl0
    ]

    lk=0
    lk=lk+ -1/2*log((a_0^2+a_1^2+a_2^2)^2-(a_0*a_1+a_1*a_2)^2)-1/2*tr(MA1inv*S1)
    lk=lk+ -1/2*log((a_0^2+a_1^2+a_2^2)^2-(a_0*a_2)^2)-1/2*tr(MA2inv*S2)
    
    

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
    fa1_expr = eval(Meta.parse(fa1_str) )
    fa2_expr = eval(Meta.parse(fa2_str))

    F = System([fa0_expr,fa1_expr,fa2_expr])

    resultcl2 = solve(F)
    
#    l0= (deterM^2/det(MA1)^2)*(tr(M1inv*SM0)*det(MA1)-tr(M1inv*SM0*M1inv*S1))+(deterM^2/det(M2)^2)*(tr(M2inv*SM0)*det(M2)-tr(M2inv*SM0*M2inv*S2))
#    l1= (tr(M1inv*SMH)*det(MA1)-tr(M1inv*SMH*M1inv*S1))
#    l2= (tr(M2inv*SMH)*det(MA2)-tr(M2inv*SMH*M2inv*S2))
 
#    f0=simplify(l0,expand=true)
#    f1=simplify(l1,expand=true)
#    f2=simplify(l2,expand=true)
 
#    @var g0 g1 g2 g3 g4 g5 g6

#    fa0 = string(f0)
#    fa1 = string(f1)
#    fa2 = string(f2)
    #fa3 = string(f3)
    #fa4 = string(f4)
    #fa5 = string(f5)
    #fa6 = string(f6)

#    fa0_str = replace(fa0, "g_0" => "g0","g_1" => "g1","g_2" => "g2")
#    fa1_str = replace(fa1, "g_0" => "g0","g_1" => "g1","g_2" => "g2")
#    fa2_str = replace(fa2, "g_0" => "g0","g_1" => "g1","g_2" => "g2")


#    fa0_expr = eval(Meta.parse(fa0_str))
#    fa1_expr = eval(Meta.parse(fa1_str))
#    fa2_expr = eval(Meta.parse(fa2_str))

#    F = System([fa0_expr,fa1_expr,fa2_expr])

#    resultcl2 = solve(F)
    
    z=solutions(resultcl2)
    tot=length(z)


    soltotales=[]
    solr=z
    soltotales=z

#    for ind in 1:tot
#        fg0=a0^2+a1^2+a2^2-solr[ind][1]
#        fg1=a0*a1+a1*a2-solr[ind][2]
#        fg2=a0*a2-solr[ind][3]


#        fag0 = string(fg0)
#        fag1 = string(fg1)
#        fag2 = string(fg2)


#        fag0_str = replace(fag0, "a_0" => "a0","a_1" => "a1","a_2" => "a2")
#        fag1_str = replace(fag1, "a_0" => "a0","a_1" => "a1","a_2" => "a2")
#        fag2_str = replace(fag2, "a_0" => "a0","a_1" => "a1","a_2" => "a2")


  #      fag0_expr = eval(Meta.parse(fag0_str))
 #       fag1_expr = eval(Meta.parse(fag1_str))
#        fag2_expr = eval(Meta.parse(fag2_str))


#        F = System([fg0,fg1,fg2])
    
#        resultx0 = solve(F)

#        for i in 1:8
 #           println(solutions(resultx0)[i])
#            push!(soltotales,solutions(resultx0)[i])
#        end
#    end

    solfinal=[]
    ns=length(soltotales)


    for i in 1:ns
       solu=soltotales[i]
        flag=true
        for j in 1:3
           if imag(solu[j])>10^(-8)
            flag=false
         end
        end
        if real(solu[1])<0
           flag=false
        end
        if abs(solu[1]) <= abs(solu[3]) 
            flag = false
        end
        if flag 
           # println(i)
            push!(solfinal,[real(solu[1]),real(solu[2]),real(solu[3])])    
        end

    end
    
    
    if length(solfinal) == 0
    for i in 1:ns
        solu = soltotales[i]
            if real(solu[1]) >= 0 && abs(solu[1]) > abs(solu[3])
                push!(solfinal, [real(solu[1]), real(solu[2]), real(solu[3])])
            end
        end
    end
    

    totsolf=length(solfinal)

   indtot=[]
    max=-Inf
    

    for i in 1:totsolf
          solt=solfinal[i]
          valores = Dict(a_0 => solt[1], a_1 => solt[2], a_2 => solt[3])
          exp_num = substitute(lk, valores)
          
          if exp_num-max>10^(-4)
            indtot=[]
            push!(indtot,i)
            max=exp_num
#            println([i,max])
          elseif abs(exp_num-max)<=10^(-10)
             push!(indtot,i)
          end

    end
    solf2=[]
    for s2 in indtot
        push!(solf2,solfinal[s2])
    end



    

    solf2 = []
    for s2 in indtot
        push!(solf2, solfinal[s2])
    end

    # If there are more than 2 solutions, keep the one closest to the target vector
    if length(solf2) > 1
        # Compute Euclidean distance to target for each solution
       distances_to_target = [sqrt(sum((sol .- [1,0.5,-0.3]).^2)) for sol in solf2]
    
        # Find the index of the solution with minimum distance
        idx_closest = argmin(distances_to_target)
    
        # Keep only that solution
        solf2 = [solf2[idx_closest]]
    end

    push!(estim, solf2)


end


estim

using Statistics   # For mean() and std()


all_solutions = vcat(estim...)  # Concatenate all arrays in `estim` into one array

solution_matrix = hcat(all_solutions...)'  # Transpose so rows = solutions, columns = parameters

mean_values = mean(solution_matrix, dims=1)         # Mean of each column
std_values = std(solution_matrix, dims=1)           # Standard deviation of each column


println("Mean of each parameter (a_0, a_1, a_2): ", mean_values)
println("Standard deviation of each parameter (a_0, a_1, a_2): ", std_values)
