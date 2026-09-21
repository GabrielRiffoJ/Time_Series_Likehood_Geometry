using HomotopyContinuation
using Symbolics
using LinearAlgebra
using DelimitedFiles
using Statistics

@variables a_1

Samples = readdlm("samples.txt", '\t', Float64, '\n', skipstart=1)
nrep = size(Samples, 1)
n = size(Samples, 2)

M = zeros(Num, n, n)

for i in 1:n
    M[i, i] = a_1^2+1
end

for i in 1:(n-1)
    M[i, i+1] = -a_1
    M[i+1, i] = -a_1
end

D1 = zeros(Num, n, n)

for i in 1:n
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

M[1,1]=1
M[n,n]=1

D1[1,1]=0
D1[n,n]=0

deter=1-(a_1^2)

Minvdet=Minv*deter

for i in 1:n
    for j in 1:n
        Minvdet[i,j]=simplify(Minv[i,j]*deter,expand=true)
    end
end

@var a1

res=[]
sigma_res=[]
for muestra in 1:nrep
    Y=Samples[muestra,:]

    lik=-log(1/(1-a_1^2))-n*log(transpose(Y)*M*Y)
    siges=1/n*transpose(Y)*M*Y

    l1= tr(D1*Minvdet)*transpose(Y)*M*Y-n*transpose(Y)*D1*Y*deter

    f1=simplify(l1,expand=true)

    fa1 = string(f1)
    fa1_str = replace(fa1,"a_1" => "a1")
    fa1_expr = eval(Meta.parse(fa1_str))

    F = System([fa1_expr])

    resultcl2 = solve(F)

    sf=real(solutions(resultcl2)[1])
    lkf=-Inf
    for i in 1:length(solutions(resultcl2))
        s=solutions(resultcl2)[i]
        if (abs(imag(s)[1])<0.001 && abs(real(s)[1])<1)
            lik_val = substitute(lik, a_1 => real(s)[1])
            if lik_val>lkf
                sf=real(s)[1]
                lkf=lik_val
            end
        end
    end
    sigf=substitute(siges,a_1=>sf)
    push!(res,sf)
    push!(sigma_res,sigf)
end

res=Float64.(res)
sigma_res=Float64.(sigma_res)

println(mean(res))
println(mean(res)-0.5)
println(sqrt(mean((res .- mean(res)).^2)))

println(mean(sigma_res))
println(mean(sigma_res)-1)
println(sqrt(mean((sigma_res .- mean(sigma_res)).^2)))
