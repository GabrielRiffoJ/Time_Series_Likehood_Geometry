using HomotopyContinuation
using Symbolics
using LinearAlgebra

@variables g_0 g_1 g_2 g_3 g_4 g_5 g_6


Pkg.add("Primes")
using Primes

prime = collect(primes(1, 1987))[1:300]

#YTime = Vector{Float64}(undef, 100)

theta1 = 0.6
theta2 = -0.3
sigma  = 1.0

n = 300


ep = rand(Normal(0, sigma), n + 2)


Ym = zeros(n)
for t in 1:n
    Ym[t] = ep[t+2] + theta1 * ep[t+1] + theta2 * ep[t]
end


n=3

YI=prime
YF=  1.0 ./ prime

#first approach
Y=YI[1:n]



q=6
M1 = zeros(Num, 2, 2)

M1=[ g_0   g_1
     g_1   g_0
]


M2=[ g_0   g_2
     g_2   g_0
]

M3=[ g_0   g_3
     g_3   g_0
]

M4=[ g_0   g_4
     g_4   g_0
]


M5=[ g_0   g_5
     g_5   g_0
]


M6=[ g_0   g_6
     g_6   g_0
]


M7=[ g_0   g_7
     g_7   g_0
]



SM0=[1   0
     0   1
]

SMH=[0   1
     1   0
]

S1=[ 
1 1/2
1/2 1
]

S2=[
1 1/3
1/3 1
]



S3=[
1 1/5
1/5  1
]

S4=[
1 1/7
1/7 1
]

S5=[
1  1/11
1/11   1
]

S6=[
1 1/13
1/13 1
]

S7=[
#1   1/17
#1/17 1

]




deterM=1
deterM=deterM*det(M1)*det(M2)*det(M3)*det(M4)*det(M5)*det(M6)*det(M7)
M1inv=inv(M1)*det(M1)
M2inv=inv(M2)*det(M2)
M3inv=inv(M3)*det(M3)
M4inv=inv(M4)*det(M4)
M5inv=inv(M5)*det(M5)
M6inv=inv(M6)*det(M6)
M6inv=inv(M6)*det(M7)



l0= (deterM^2/det(M1)^2)*(tr(M1inv*SM0)*det(M1)-tr(M1inv*SM0*M1inv*S1))+(deterM^2/det(M2)^2)*(tr(M2inv*SM0)*det(M2)-tr(M2inv*SM0*M2inv*S2))+(deterM^2/det(M3)^2)*(tr(M3inv*SM0)*det(M3)-tr(M3inv*SM0*M3inv*S3))+(deterM^2/det(M4)^2)*(tr(M4inv*SM0)*det(M4)-tr(M4inv*SM0*M4inv*S4))+(deterM^2/det(M5)^2)*(tr(M5inv*SM0)*det(M5)-tr(M5inv*SM0*M5inv*S5))+(deterM^2/det(M6)^2)*(tr(M6inv*SM0)*det(M6)-tr(M6inv*SM0*M6inv*S6))+(deterM^2/det(M7)^2)*(tr(M7inv*SM0)*det(M7)-tr(M7inv*SM0*M7inv*S7))
l1= (tr(M1inv*SMH)*det(M1)-tr(M1inv*SMH*M1inv*S1))
l2= (tr(M2inv*SMH)*det(M2)-tr(M2inv*SMH*M2inv*S2))
l3= (tr(M3inv*SMH)*det(M3)-tr(M3inv*SMH*M3inv*S3))
l4= (tr(M4inv*SMH)*det(M4)-tr(M4inv*SMH*M4inv*S4))
l5= (tr(M5inv*SMH)*det(M5)-tr(M5inv*SMH*M5inv*S5))
l6= (tr(M6inv*SMH)*det(M6)-tr(M6inv*SMH*M6inv*S6))
l7= (tr(M7inv*SMH)*det(M7)-tr(M7inv*SMH*M7inv*S7))


f0=simplify(l0,expand=true)
f1=simplify(l1,expand=true)
f2=simplify(l2,expand=true)
f3=simplify(l3,expand=true)
f4=simplify(l4,expand=true)
f5=simplify(l5,expand=true)
f6=simplify(l6,expand=true)
f7=simplify(l6,expand=true)


@var g0 g1 g2 g3 g4 g5 g6 g7

fa0 = string(f0)
fa1 = string(f1)
fa2 = string(f2)
fa3 = string(f3)
fa4 = string(f4)
fa5 = string(f5)
fa6 = string(f6)
fa7 = string(f7)



fa0_str = replace(fa0, "g_0" => "g0","g_1" => "g1","g_2" => "g2","g_3" => "g3","g_4" => "g4","g_5" => "g5","g_6" => "g6","g_7" => "g7")
fa1_str = replace(fa1, "g_0" => "g0","g_1" => "g1","g_2" => "g2","g_3" => "g3","g_4" => "g4","g_5" => "g5","g_6" => "g6","g_7" => "g7")
fa2_str = replace(fa2, "g_0" => "g0","g_1" => "g1","g_2" => "g2","g_3" => "g3","g_4" => "g4","g_5" => "g5","g_6" => "g6","g_7" => "g7")
fa3_str = replace(fa3, "g_0" => "g0","g_1" => "g1","g_2" => "g2","g_3" => "g3","g_4" => "g4","g_5" => "g5","g_6" => "g6","g_7" => "g7")
fa4_str = replace(fa4, "g_0" => "g0","g_1" => "g1","g_2" => "g2","g_3" => "g3","g_4" => "g4","g_5" => "g5","g_6" => "g6","g_7" => "g7")
fa5_str = replace(fa5, "g_0" => "g0","g_1" => "g1","g_2" => "g2","g_3" => "g3","g_4" => "g4","g_5" => "g5","g_6" => "g6","g_7" => "g7")
fa6_str = replace(fa6, "g_0" => "g0","g_1" => "g1","g_2" => "g2","g_3" => "g3","g_4" => "g4","g_5" => "g5","g_6" => "g6","g_7" => "g7")


fa0_expr = eval(Meta.parse(fa0_str))
fa1_expr = eval(Meta.parse(fa1_str))
fa2_expr = eval(Meta.parse(fa2_str))
fa3_expr = eval(Meta.parse(fa3_str))
fa4_expr = eval(Meta.parse(fa4_str))
fa5_expr = eval(Meta.parse(fa5_str))
fa6_expr = eval(Meta.parse(fa6_str))
fa7_expr = eval(Meta.parse(fa7_str))



F = System([fa0_expr,fa1_expr,fa2_expr,fa3_expr,fa4_expr,fa5_expr,fa6_expr,fa7_expr])

#F = System([fa0_expr,fa1_expr,fa2_expr,fa3_expr,fa4_expr,fa5_expr])



t1 = Sys.time()

resultcl2 = solve(F)

t2=Sys.time()

t2-t1



#S1+S2+S3+S4+S5+S6

#TS=[
 #   1   1/2 1/3 1/5 1/7 1/11 
  #  1/2 1   1/2 1/3 1/5 1/7
#    1/3 1/2 1 1/2 1/3 1/5 
#    1/5 1/3 1/2 1 1/2 1/3
#    1/7 1/5 1/3 1/2 1 1/2
#    1/11 1/7 1/5 1/3 1/2 1   
#]
