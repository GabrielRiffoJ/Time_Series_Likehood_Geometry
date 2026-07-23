using HomotopyContinuation
using Symbolics
using LinearAlgebra

#@variables a_1 a_2 a_3 a_4 a_5 a_6 a_7 a_8 a_9 a_10

@variables a_1

n=12

M = zeros(Num, n, n)


for i in 1:n
    M[i, i] = a_1^2+1
end

for i in 1:(n-1)
    M[i, i+1] = -a_1
    M[i+1, i] = -a_1
end


D1 = zeros(Num, n, n)

for i in 1:(n)
    D1[i, i] = 2*a_1
end

for i in 1:(n-1)
    D1[i, i+1] = -1
    D1[i+1, i] = -1
end




Minv = zeros(Num, n, n)
    

for h in 0:(n-1)
    for i in 1:(n-h)
        Minv[i+h, i] = a_1^(h)/(1-a_1^2)
        Minv[i, i+h] = a_1^(h)/(1-a_1^2)
    end
end

Minv

#YI=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73]
#YF=[1/3,1/5,1/7,1/11,1/13,1/17,1/19,1/23,1/29,1/31,1/37]
#Y=YI[1:n]

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



#Minv=inv(M)
deter=1-(a_1^2)


Minvdet=Minv*deter


for i in 1:n
    for j in 1:n
        Minvdet[i,j]=simplify(Minv[i,j]*deter,expand=true)
    end
end



@variables x y
@variables a_1 a_2
@var a1


res=[]
sigma_res=[]
for muestra in 1:100
    Y=Samples[muestra]

    lik=-log(1/(1-a_1^2))-n*log(transpose(Y)*M*Y)
    siges=1/n*transpose(Y)*M*Y

    l1= tr(D1*Minvdet)*transpose(Y)*M*Y-n*transpose(Y)*D1*Y*deter

    f1=simplify(l1,expand=true)

    fa1 = string(f1)

    @var a0 a1 a2

    fa1_str = replace(fa1,"a_1" => "a1")

    fa1_expr = eval(Meta.parse(fa1_str))

    F = System([fa1_expr])

    resultcl2 = solve(F)
    
    sf=real(solutions(resultcl2)[1])
    lkf=-Inf
    for i in 1:length(solutions(resultcl2))
        s=solutions(resultcl2)[i]
        if (abs(imag(s)[1])<0.001 &&  abs(real(s)[1])<1)
            lik_val = substitute(lik, a_1 => real(s)[1])
            if lik_val>lkf
                if abs(real(s)[1])<1
                    sf=real(s)[1] 
                    lkf=lik_val
                end
            end 
        end
    end
    sigf=substitute(siges,a_1=>sf)
    push!(res,sf)
    push!(sigma_res,sigf)
end



suma=0
for i in 1:100
    suma=suma+res[i]
end


suma_sig=0
for i in 1:100
    suma_sig=suma_sig+sigma_res[i]
end

sdt1=0
for i in 1:100
    sdt1=sdt1+(res[i]-suma/100)^2
end
sqrt(sdt1/100)

sdsig=0

for i in 1:100
    sdsig=sdsig+(sigma_res[i]-suma_sig/100)^2
end

sqrt(sdsig/100)

println(suma/100)
println(suma/100-0.5)
println(sqrt(sdt1/100))


println(suma_sig/100)
println(suma_sig/100-1)
println(sqrt(sdsig/100))

