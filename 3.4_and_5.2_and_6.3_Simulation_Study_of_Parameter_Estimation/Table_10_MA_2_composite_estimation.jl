using HomotopyContinuation
using Symbolics
using LinearAlgebra
using DelimitedFiles
using Statistics

@variables a_0 a_1 a_2
@var a0 a1 a2

# Empirical autocovariances (gamma0, gamma1, gamma2), one row per replication
G = readdlm("samples_MA2_composite.txt", '\t', Float64, '\n', skipstart=1)
nrep = size(G, 1)

# Covariance matrices of the pairs (Y_t, Y_{t+1}) and (Y_t, Y_{t+2})
M1 = [ a_0^2+a_1^2+a_2^2   a_0*a_1+a_1*a_2
       a_0*a_1+a_1*a_2     a_0^2+a_1^2+a_2^2 ]

M2 = [ a_0^2+a_1^2+a_2^2   a_0*a_2
       a_0*a_2             a_0^2+a_1^2+a_2^2 ]

# Derivatives of M1 with respect to a_0, a_1, a_2
a0M1 = [ 2*a_0  a_1
         a_1    2*a_0 ]
a1M1 = [ 2*a_1      a_0+a_2
         a_0+a_2    2*a_1 ]
a2M1 = [ 2*a_2  a_1
         a_1    2*a_2 ]

# Derivatives of M2 with respect to a_0, a_1, a_2
a0M2 = [ 2*a_0  a_2
         a_2    2*a_0 ]
a1M2 = [ 2*a_1  0
         0      2*a_1 ]
a2M2 = [ 2*a_2  a_0
         a_0    2*a_2 ]

detM1 = det(M1)
detM2 = det(M2)

# Adjugates of the 2x2 matrices (inverse = adjugate / determinant)
M1inv = [ M1[2,2] -M1[1,2]
          -M1[2,1] M1[1,1] ]
M2inv = [ M2[2,2] -M2[1,2]
          -M2[2,1] M2[1,1] ]

estim = Vector{Vector{Float64}}()
fallidas = 0

for act in 1:nrep

    println(act)

    gcl0 = G[act, 1]
    gcl1 = G[act, 2]
    gcl2 = G[act, 3]

    # Empirical matrices for lags 1 and 2
    S1 = [ gcl0 gcl1
           gcl1 gcl0 ]
    S2 = [ gcl0 gcl2
           gcl2 gcl0 ]

    # Composite log-likelihood (up to constants)
    lk = -1/2*log(detM1) - 1/2*tr((M1inv/detM1)*S1) - 1/2*log(detM2) - 1/2*tr((M2inv/detM2)*S2)

    # Score of each lag term multiplied by det^2 to clear denominators
    f1a0 = detM1*(-0.5*tr(a0M1*M1inv)) - 0.5*tr((-(M1inv)*a0M1*M1inv)*S1)
    f1a1 = detM1*(-0.5*tr(a1M1*M1inv)) - 0.5*tr((-(M1inv)*a1M1*M1inv)*S1)
    f1a2 = detM1*(-0.5*tr(a2M1*M1inv)) - 0.5*tr((-(M1inv)*a2M1*M1inv)*S1)

    f2a0 = detM2*(-0.5*tr(a0M2*M2inv)) - 0.5*tr((-(M2inv)*a0M2*M2inv)*S2)
    f2a1 = detM2*(-0.5*tr(a1M2*M2inv)) - 0.5*tr((-(M2inv)*a1M2*M2inv)*S2)
    f2a2 = detM2*(-0.5*tr(a2M2*M2inv)) - 0.5*tr((-(M2inv)*a2M2*M2inv)*S2)

    # Total score as a polynomial system (common denominator det1^2 * det2^2)
    fa0ex = Symbolics.simplify(detM2^2*f1a0 + detM1^2*f2a0, expand=true)
    fa1ex = Symbolics.simplify(detM2^2*f1a1 + detM1^2*f2a1, expand=true)
    fa2ex = Symbolics.simplify(detM2^2*f1a2 + detM1^2*f2a2, expand=true)

    fa0 = string(fa0ex)
    fa1 = string(fa1ex)
    fa2 = string(fa2ex)

    # Rename Symbolics variables to HomotopyContinuation variables
    fa0_str = replace(fa0, "a_0" => "a0", "a_1" => "a1", "a_2" => "a2")
    fa1_str = replace(fa1, "a_0" => "a0", "a_1" => "a1", "a_2" => "a2")
    fa2_str = replace(fa2, "a_0" => "a0", "a_1" => "a1", "a_2" => "a2")

    fa0_expr = eval(Meta.parse(fa0_str))
    fa1_expr = eval(Meta.parse(fa1_str))
    fa2_expr = eval(Meta.parse(fa2_str))

    F = System([fa0_expr, fa1_expr, fa2_expr])

    # Track all complex solutions of the score equations
    resultcl2 = solve(F)

    soltotales = solutions(resultcl2)

    solfinal = Vector{Vector{Float64}}()

    # Keep real solutions with a_0 >= 0 (removes the global sign symmetry)
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

    # Keep only the invertible representation: both roots of
    # a_0 z^2 + a_1 z + a_2 lie inside the unit disk
    # (equivalently, both roots of a_0 + a_1 x + a_2 x^2 lie outside it).
    # This rule does not use the true parameter values.
    solfinal = filter(solfinal) do s
        d = sqrt(Complex(s[2]^2 - 4*s[1]*s[3]))
        z1 = (-s[2] + d)/(2*s[1])
        z2 = (-s[2] - d)/(2*s[1])
        abs(z1) < 1 && abs(z2) < 1
    end

    # Among the remaining critical points, keep the one with the largest likelihood
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
        # Replications with no invertible real solution are counted, not hidden
        global fallidas += 1
        println("no invertible solution for replication: ", act)
    end
end

# Summary statistics over the replications that produced an estimate
E = reduce(hcat, estim)'
verdad = [1.0, 0.5, -0.3]
medias = vec(mean(E, dims=1))

println(medias)
println(medias .- verdad)
println([sqrt(mean((E[:, j] .- medias[j]).^2)) for j in 1:3])
println(fallidas)
