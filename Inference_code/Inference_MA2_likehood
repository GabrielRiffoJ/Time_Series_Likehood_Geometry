using HomotopyContinuation
using Symbolics
using LinearAlgebra

@variables a_0 a_1 a_2



n=10

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


detM=det(M)
detMex=simplify(detM,expand=true)

Minv2=(detM*Minv)
Minv2ex=Minv2


Ma = zeros(Num, n, n)
    

for i in 1:n
    Ma[i, i] = a_0^2+a_1^2+a_2^2
end
    

for i in 2:n
    Ma[i, i-1] = a_0*a_1+a_1*a_2
    Ma[i-1, i] = a_0*a_1+a_1*a_2
end
    

for i in 3:n
    Ma[i, i-2] = a_0*a_2
    Ma[i-2, i] = a_0*a_2
end


#@variables y1 y2 y3 y4 y5 y6 y7 y8 y9 y10
#Y2=[y1, y2, y3, y4, y5, y6, y7, y8, y9, y10]



Mainv=inv(Ma)
detera=det(Ma)


deter=det(M)

Minv=inv(M)


Minvdet=Minv*deter

@variables a0 a1 a2

t1=Sys.time()

estim=[]

nmuestras=2

for mact in 1:nmuestras

    println(mact)
    Y=muestras[mact,1:10]
    lk= -1/2*log(detera)-1/2*transpose(Y)*Mainv*Y
    
    l0= tr(Minvdet*Ma0)*deter-transpose(Y)*Minvdet*Ma0*Minvdet*Y
    l1= tr(Minvdet*Ma1)*deter-transpose(Y)*Minvdet*Ma1*Minvdet*Y
    l2= tr(Minvdet*Ma2)*deter-transpose(Y)*Minvdet*Ma2*Minvdet*Y


    f0=simplify(l0,expand=true)
    f1=simplify(l1,expand=true)
    f2=simplify(l2,expand=true)


    fa0 = string(f0)
    fa1 = string(f1)
    fa2 = string(f2)


    fa0_str = replace(fa0, "a_0" => "a0","a_1" => "a1","a_2" => "a2")
    fa1_str = replace(fa1, "a_0" => "a0","a_1" => "a1","a_2" => "a2")
    fa2_str = replace(fa2, "a_0" => "a0","a_1" => "a1","a_2" => "a2")


    fa0_expr = eval(Meta.parse(fa0_str))
    fa1_expr = eval(Meta.parse(fa1_str))
    fa2_expr = eval(Meta.parse(fa2_str))

    F = System([fa0_expr,fa1_expr,fa2_expr])

#    t1 = time()

    resultcl2 = solve(F)


    #t2=Sys.time()

 #   t2-t1


#t3=Sys.time()
    z=solutions(resultcl2)
    tot=length(z)
    

    soltotales=z


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
        if flag 
#            println(i)
            push!(solfinal,[real(solu[1]),real(solu[2]),real(solu[3])])    
        end

    end

    totsolf=length(solfinal)

    indtot=0
    max=-Inf

    for i in 1:totsolf
          solt=solfinal[i]
          valores = Dict(a_0 => solt[1], a_1 => solt[2], a_2 => solt[3])
          exp_num = substitute(lk, valores)

          if exp_num>max
            indtot=i
            max=exp_num
#            println([i,max])
          end
    end
    push!(estim,solfinal[indtot])
end


ind
indtot

abs(so12[1])+abs(so12[2])+abs(so12[3])




