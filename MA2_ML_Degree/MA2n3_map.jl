using HomotopyContinuation
using Symbolics
using LinearAlgebra

@variables a_0 a_1 a_2


n=3

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

using Random, Primes


primes_list = primes(2, 2000)  
primes_list = primes_list[1:300]


reciprocals = 1.0 ./ primes_list


m = 300
e = randn(m)
theta1, theta2 = 0.5, 0.3
ma2 = e .+ theta1 .* [0.0; e[1:end-1]] .+ theta2 .* [0.0; 0.0; e[1:end-2]]


primes_list, reciprocals, ma2


Y=primes_list[1:n]

deter=det(M)

Minv=inv(M)

Minv*M

Minvdet=Minv*deter




# ==========================================
# Auxiliar Function
# ==========================================

# Discriminant: y1^4 - 4y1^2*y2^2 + y2^4 + 6y1*y2^2*y3 - y1^2*y3^2 - 4y2^2*y3^2 + y3^4
sqrt_term(y1, y2, y3) = sqrt(
    y1^4 - 4*y1^2*y2^2 + y2^4 +
    6*y1*y2^2*y3 - y1^2*y3^2 -
    4*y2^2*y3^2 + y3^4
)

# Denominator: (y1 - 2y2 + y3)*(y1 + 2y2 + y3)
den(y1, y2, y3) =
    (y1 - 2*y2 + y3) *
    (y1 + 2*y2 + y3)

# ==========================================
# γ0
# ==========================================
function gamma_0(y1, y2, y3; branch=:plus)
    # ----------------------------------
    # Special case I:  y1 + 2y2 + y3 = 0
    # ----------------------------------
    if abs(y1 + 2*y2 + y3) < 1e-12
        g0 = (y3^2 + 2*y2*y3 - 5*y2^2) / 6
        return g0
    end
    # ----------------------------------
    # Special case II: y1 - 2y2 + y3 = 0
    # ----------------------------------
    if abs(y1 - 2*y2 + y3) < 1e-12
        g0 = (y3^2 - 2*y2*y3 - 5*y2^2) / 6
        return g0
    end
    # ----------------------------------
    # General case
    # ----------------------------------
    s = sqrt_term(y1, y2, y3)
    d = den(y1, y2, y3)

    base =
        y1^4 - 10*y2^4 +
        2*y1^3*y3 +
        2*y1*y2^2*y3 +
        2*y1*y2*y3^3 +       
        2*y1*y3^3 +
        y3^4


    corr = ((y1 + y3)^2 + 2*y2^2) * s

    if branch == :plus
        return (base + corr) / (6*d)
    elseif branch == :minus
        return (base - corr) / (6*d)
    else
        error("branch must be :plus or :minus")
    end
end

# ==========================================
# γ1
# ==========================================
function gamma_1(y1, y2, y3; branch=:plus)
    # ----------------------------------
    # Special case I:  y1 + 2y2 + y3 = 0
    # ----------------------------------
    if abs(y1 + 2*y2 + y3) < 1e-12
        return y2^2
    end
    # ----------------------------------
    # Special case II: y1 - 2y2 + y3 = 0
    # ----------------------------------
    if abs(y1 - 2*y2 + y3) < 1e-12
        return -y2^2
    end
    # ----------------------------------
    # General case
    # ----------------------------------
    s = sqrt_term(y1, y2, y3)
    d = den(y1, y2, y3)

    core = y1^2 - 3*y2^2 + y1*y3 + y3^2

    if branch == :plus
        num = y2*(y1 + y3)*(core - s)
    elseif branch == :minus
        num = y2*(y1 + y3)*(core + s)
    else
        error("branch must be :plus or :minus")
    end

    # Denominator: 2*(y1-2y2+y3)*(y1+2y2+y3)  =  2*d
    return num / (2*d)
end

# ==========================================
# γ2
# ==========================================
function gamma_2(y1, y2, y3; branch=:plus)
    # ----------------------------------
    # Special case I:  y1 + 2y2 + y3 = 0
    # ----------------------------------
    if abs(y1 + 2*y2 + y3) < 1e-12
        y1_val = -2*y2 - y3      # y1 = -2y2 - y3
        return -(y3^2 + 2*y2*y3 + 3*y2^2) / 2
    end
    # ----------------------------------
    # Special case II: y1 - 2y2 + y3 = 0
    # ----------------------------------
    if abs(y1 - 2*y2 + y3) < 1e-12
        return -(y3^2 - 2*y2*y3 + 3*y2^2) / 2
    end
    # ----------------------------------
    # General case
    # ----------------------------------
    s = sqrt_term(y1, y2, y3)
    d = den(y1, y2, y3)

    base =
        y1^2*y2^2 -
        2*y2^4 +
        y1^3*y3 -
        4*y1*y2^2*y3 +
        2*y1^2*y3^2 +
        y2^2*y3^2 +
        y1*y3^3


    # Denominator: 2*(y1-2y2+y3)*(y1+2y2+y3)  =  2*d
    if branch == :plus
        return (base - 2*y2*s) / (2*d)
    elseif branch == :minus
        return (base + 2*y2*s) / (2*d)
    else
        error("branch must be :plus or :minus")
    end
end



@variables y_[1:10]


#Y_vec = [Y[i] for i in 1:n]
#Y_vec=Y[1:n]
#Y_vec=[y_[1],y_[2],y_[3],y_[4],y_[5],y_[6],y_[7],y_[8],y_[9],y_[10]]

#Y=Y_vec[1:n]

y=[1,0.5,1]
s1a0=gamma_0(y[1],y[2],y[3], branch=:plus)
s1a1=gamma_1(y[1],y[2],y[3], branch=:plus)
s1a2=gamma_2(y[1],y[2],y[3], branch=:plus)

s2a0=gamma_0(y[1],y[2],y[3], branch=:minus)
s2a1=gamma_1(y[1],y[2],y[3], branch=:minus)
s2a2=gamma_2(y[1],y[2],y[3], branch=:minus)



subs1 = Dict(a_0 => s1a0, a_1 => s1a1, a_2 => s1a2)
subs2 = Dict(a_0 => s2a0, a_1 => s2a1, a_2 => s2a2)

Minv_s1 = copy(Minv)
Minv_s2 = copy(Minv)

for i in 1:n
    for j in 1:n
        Minv_s1[i, j] = Symbolics.substitute(Minv[i, j], subs1)
        Minv_s2[i, j] = Symbolics.substitute(Minv[i, j], subs2)
    end
end

deter_s1=Symbolics.substitute(deter, subs1)
deter_s2=Symbolics.substitute(deter, subs2)

sol=0

if deter_s1 < 0 && deter_s2 < 0
    error("Determinant is non-positive, cannot compute likelihood.")
elseif deter_s1<0
    sol=2
elseif deter_s2<0
    sol=1
else
    f_like=sqrt(deter_s1)^(-1)*exp(-0.5* transpose(y)*Minv_s1*y)
    f_like2=sqrt(deter_s2)^(-1)*exp(-0.5* transpose(y)*Minv_s2*y)
    if f_like > f_like2
        sol=1
    elseif f_like2 > f_like
        sol=2
    else
        sol=0
    end
end


f_like=sqrt(deter_s1)^(-1)*exp(-0.5* transpose(y)*Minv_s1*y)
f_like2=sqrt(deter_s2)^(-1)*exp(-0.5* transpose(y)*Minv_s2*y)








using Symbolics
using LinearAlgebra
using Plots


function classify(y1,y2,y3)


    if abs(y1 - 2y2 + y3) < 1e-12 ||    abs(y1 + 2y2 + y3) < 1e-12

            return -1
    end

    if y1^4 - 4*y1^2*y2^2 + y2^4 +   6*y1*y2^2*y3 - y1^2*y3^2 - 4*y2^2*y3^2 + y3^4 < 0
            return 0
    end

    y = [y1,y2,y3]

    # -------------------------
    # Solution 1
    # -------------------------

    s1a0 = gamma_0(y1,y2,y3, branch=:plus)
    s1a1 = gamma_1(y1,y2,y3, branch=:plus)
    s1a2 = gamma_2(y1,y2,y3, branch=:plus)

    # -------------------------
    # Solution 2
    # -------------------------

    s2a0 = gamma_0(y1,y2,y3, branch=:minus)
    s2a1 = gamma_1(y1,y2,y3, branch=:minus)
    s2a2 = gamma_2(y1,y2,y3, branch=:minus)

    # -------------------------
    # Substitutions
    # -------------------------

    subs1 = Dict(a_0 => s1a0, a_1 => s1a1, a_2 => s1a2)
    subs2 = Dict(a_0 => s2a0, a_1 => s2a1, a_2 => s2a2)

    # -------------------------
    # Evaluate determinants
    # -------------------------

    deter_s1 = Symbolics.value(Symbolics.substitute(deter, subs1))
    deter_s2 = Symbolics.value(Symbolics.substitute(deter, subs2))

    # -------------------------
    # Invalid cases
    # -------------------------

    if deter_s1 < 0 && deter_s2 < 0
        return 0
    elseif deter_s1 < 0
        return 2
    elseif deter_s2 < 0
        return 1
    end

    # -------------------------
    # Build matrices
    # -------------------------

    Minv_s1 = Matrix{Float64}(undef, n, n)
    Minv_s2 = Matrix{Float64}(undef, n, n)

    for i in 1:n
        for j in 1:n

            Minv_s1[i,j] =
                Symbolics.value(
                    Symbolics.substitute(Minv[i,j], subs1)
                )

            Minv_s2[i,j] =
                Symbolics.value(
                    Symbolics.substitute(Minv[i,j], subs2)
                )
        end
    end

    # -------------------------
    # Likelihoods
    # -------------------------

    f_like =
        deter_s1^(-0.5) *
        exp(-0.5 * transpose(y) * Minv_s1 * y)

    f_like2 =
        deter_s2^(-0.5) *
        exp(-0.5 * transpose(y) * Minv_s2 * y)

    # -------------------------
    # Final classification
    # -------------------------

    if f_like > f_like2
        return 1

    elseif f_like2 > f_like
        return 2

    else
        return 0
    end
end




xvals = range(-50,50,length=200)
yvals = range(-50,50,length=200)

z_fixed = 2

Z = [
    classify(x,y,z_fixed)
    for y in yvals, x in xvals
]

# ==========================================
# GRAFIC
# ==========================================

Zplot = Float64.(Z)

VY3=heatmap(
    xvals,
    yvals,
    Zplot,

    color = cgrad(
        [:green,:black, :red, :blue],
        4,
        categorical = true
    ),

    clims = (-1, 2),

    colorbar = false,
    xlabel = "y₁",
    ylabel = "y₂",
    title = ""
)



plot!(
    xvals,
    -0.5*(-xvals .- z_fixed),
    color = :green,
    lw = 4,
    label = false,
    xlims = (first(xvals), last(xvals)),
    ylims = (first(yvals), last(yvals))
)

plot!(
    xvals,
    0.5*(-xvals .- z_fixed),
    color = :green,
    lw = 4,
    label = false,
    xlims = (first(xvals), last(xvals)),
    ylims = (first(yvals), last(yvals))
)


VY3
savefig(VY3, "MAn3_Y_3_vf.png")



zvals = range(-50,50,length=200)
yvals = range(-50,50,length=200)

x_fixed = 2

Z = [
    classify(x_fixed,y,z)
    for z in zvals, y in yvals
]

# ==========================================
# GRAFIC
# ==========================================

Zplot = Float64.(Z)

VY1=heatmap(
    yvals,
    zvals,
    Zplot,

    color = cgrad(
        [:green, :black, :red, :blue],
        4,
        categorical = true
    ),

    clims = (-1, 2),

    colorbar = false,
    xlabel = "y₂",
    ylabel = "y₃",
    title = ""
)

plot!(
    yvals,
    -2*yvals .- x_fixed,
    color = :green,
    lw = 4,
    label = false,
    xlims = (first(yvals), last(yvals)),
    ylims = (first(zvals), last(zvals))
)
plot!(
    yvals,
    2*yvals .- x_fixed,
    color = :green,
    lw = 4,
    label = false,
    xlims = (first(yvals), last(yvals)),
    ylims = (first(zvals), last(zvals))
)

VY1


savefig(VY1, "MAn3_Y_1_vf.png")



xvals = range(-50,50,length=200)
yvals = range(-50,50,length=200)

y_fixed = 2

Z = [
    classify(x,y_fixed,z)
    for x in xvals, z in yvals
]

# ==========================================
# GRAFIC
# ==========================================

Zplot = Float64.(Z)

VY2=heatmap(
    xvals,
    zvals,
    Zplot,

    color = cgrad(
        [:green,:black, :red, :blue],
        4,
        categorical = true
    ),

    clims = (-1, 2),

    colorbar = false,
    xlabel = "y₁",
    ylabel = "y₃",
    title = ""
)

plot!(
    xvals,
    -xvals .+ 4,
    color = :green,
    lw = 4,
    label = false,
    xlims = (first(xvals), last(xvals)),
    ylims = (first(zvals), last(zvals))
)
plot!(
    xvals,
    -xvals .- 4,
    color = :green,
    lw = 4,
    label = false,
    xlims = (first(xvals), last(xvals)),
    ylims = (first(zvals), last(zvals))
)



VY2
savefig(VY2, "MAn3_Y_2_vf.png")

