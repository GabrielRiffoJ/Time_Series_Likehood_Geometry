using HomotopyContinuation
using Symbolics
using LinearAlgebra
using DelimitedFiles
using Statistics

@variables a_0 a_1 a_2

Samples = readdlm("samples_MA2.txt", '\t', Float64, '\n', skipstart=1)
nmuestras = size(Samples, 1)
n = size(Samples, 2)

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

Ma0 = Matrix{Num}(undef, n, n)
fill!(Ma0, 0)
for i in 1:n
    Ma0[i, i] = 2*a_0
    if i < n
        Ma0[i, i+1] = a_1
        Ma0[i+1, i] = a_1
    end
    if i < n-1
        Ma0[i, i+2] = a_2
        Ma0[i+2, i] = a_2
    end
end

Ma1 = Matrix{Num}(undef, n, n)
fill!(Ma1, 0)
for i in 1:n
    Ma1[i, i] = 2*a_1
    if i < n
        Ma1[i, i+1] = a_0 + a_2
        Ma1[i+1, i] = a_0 + a_2
    end
end

Ma2 = Matrix{Num}(undef, n, n)
fill!(Ma2, 0)
for i in 1:n
    Ma2[i, i] = 2*a_2
    if i < n
        Ma2[i, i+1] = a_1
        Ma2[i+1, i] = a_1
    end
    if i < n-1
        Ma2[i, i+2] = a_0
        Ma2[i+2, i] = a_0
    end
end

Minv = inv(M)
deter = det(M)
Minvdet = Minv*deter

@variables a0 a1 a2

estim = Vector{Vector{Float64}}()
fallidas = 0

for mact in 1:nmuestras

    println(mact)
    Y = Samples[mact, :]
    lk = -1/2*log(deter) - 1/2*transpose(Y)*Minv*Y

    l0 = tr(Minvdet*Ma0)*deter - transpose(Y)*Minvdet*Ma0*Minvdet*Y
    l1 = tr(Minvdet*Ma1)*deter - transpose(Y)*Minvdet*Ma1*Minvdet*Y
    l2 = tr(Minvdet*Ma2)*deter - transpose(Y)*Minvdet*Ma2*Minvdet*Y

    f0 = simplify(l0, expand=true)
    f1 = simplify(l1, expand=true)
    f2 = simplify(l2, expand=true)

    fa0 = string(f0)
    fa1 = string(f1)
    fa2 = string(f2)

    fa0_str = replace(fa0, "a_0" => "a0", "a_1" => "a1", "a_2" => "a2")
    fa1_str = replace(fa1, "a_0" => "a0", "a_1" => "a1", "a_2" => "a2")
    fa2_str = replace(fa2, "a_0" => "a0", "a_1" => "a1", "a_2" => "a2")

    fa0_expr = eval(Meta.parse(fa0_str))
    fa1_expr = eval(Meta.parse(fa1_str))
    fa2_expr = eval(Meta.parse(fa2_str))

    F = System([fa0_expr, fa1_expr, fa2_expr])

    resultcl2 = solve(F)

    soltotales = solutions(resultcl2)

    solfinal = Vector{Vector{Float64}}()

    for i in 1:length(soltotales)
        solu = soltotales[i]
        flag = true
        for j in 1:3
            if abs(imag(solu[j])) > 10^(-8)
                flag = false
            end
        end
        if real(solu[1]) < 0
            flag = false
        end
        if flag
            push!(solfinal, [real(solu[1]), real(solu[2]), real(solu[3])])
        end
    end

    solfinal = filter(solfinal) do s
        d = sqrt(Complex(s[2]^2 - 4*s[1]*s[3]))
        z1 = (-s[2] + d)/(2*s[1])
        z2 = (-s[2] - d)/(2*s[1])
        abs(z1) < 1 && abs(z2) < 1
    end

    indtot = 0
    lkmax = -Inf

    for i in 1:length(solfinal)
        solt = solfinal[i]
        valores = Dict(a_0 => solt[1], a_1 => solt[2], a_2 => solt[3])
        exp_num = substitute(lk, valores)
        if exp_num > lkmax
            indtot = i
            lkmax = exp_num
        end
    end

    if indtot > 0
        push!(estim, solfinal[indtot])
    else
        fallidas += 1
        println("muestra sin solucion invertible: ", mact)
    end
end

E = reduce(hcat, estim)'
verdad = [1.0, 0.5, -0.3]
medias = vec(mean(E, dims=1))

println(medias)
println(medias .- verdad)
println([sqrt(mean((E[:, j] .- medias[j]).^2)) for j in 1:3])
println(fallidas)




