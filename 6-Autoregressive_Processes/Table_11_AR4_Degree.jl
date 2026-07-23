using HomotopyContinuation
using Symbolics
using LinearAlgebra


@variables a_1 a_2 a_3 a_4 a_5 a_6 a_7 a_8 a_9 a_10

n=10

M = zeros(Num, n, n)
    

for i in 1:n
    M[i, i] =a_4^2+a_3^2+a_2^2+a_1^2+1
end

for i in 1:(n-1)
    M[i, i+1] = -a_1+a_1*a_2+a_2*a_3+a_3*a_4
    M[i+1, i] = -a_1+a_1*a_2+a_2*a_3+a_3*a_4
end

for i in 1:(n-2)
    M[i, i+2] = -a_2+a_1*a_3+a_2*a_4
    M[i+2, i] = -a_2+a_1*a_3+a_2*a_4
end

for i in 1:(n-3)
    M[i, i+3] = -a_3+a_1*a_4
    M[i+3, i] = -a_3+a_1*a_4
end

for i in 1:(n-4)
    M[i, i+4] = -a_4
    M[i+4, i] = -a_4
end


D1 = zeros(Num, n, n)

for i in 1:(n)
    D1[i, i] = 2*a_1
end


for i in 1:(n-1)
    D1[i, i+1] = -1+a_2
    D1[i+1, i] = -1+a_2
end

for i in 1:(n-2)
    D1[i, i+2] = a_3
    D1[i+2, i] = a_3
end

for i in 1:(n-3)
    D1[i, i+3] = a_4
    D1[i+3, i] = a_4
end

D2 = zeros(Num, n, n)

for i in 1:(n)
    D2[i, i] = 2*a_2
end


for i in 1:(n-1)
    D2[i, i+1] = a_1+a_3
    D2[i+1, i] = a_1+a_3
end

for i in 1:(n-2)
    D2[i, i+2] = -1+a_4
    D2[i+2, i] = -1+a_4
end


D3 = zeros(Num, n, n)

for i in 1:(n)
    D3[i, i] = 2*a_3
end


for i in 1:(n-1)
    D3[i, i+1] = a_2+a_4
    D3[i+1, i] = a_2+a_4
end

for i in 1:(n-2)
    D3[i, i+2] = a_1
    D3[i+2, i] = a_1
end

for i in 1:(n-3)
    D3[i, i+3] = -1
    D3[i+3, i] = -1
end

D4 = zeros(Num, n, n)

for i in 1:(n)
    D4[i, i] = 2*a_4
end


for i in 1:(n-1)
    D4[i, i+1] = a_3
    D4[i+1, i] = a_3
end

for i in 1:(n-2)
    D4[i, i+2] = a_2
    D4[i+2, i] = a_2
end

for i in 1:(n-3)
    D4[i, i+3] = a_1
    D4[i+3, i] = a_1
end

for i in 1:(n-4)
    D4[i, i+4] = -1
    D4[i+4, i] = -1
end

using Random
using Primes


nprimes = Primes.primes(2000)      
primes_300 = nprimes[1:300]             

inv_primes = 1.0 ./ primes_300


m = 300
phi = [0.4, 0.3, 0.2, 0.1]
sigma = 1.0


ar4 = zeros(m)
epsilon = sigma .* randn(m)

ar4[1] = epsilon[1]
ar4[2] = phi[1]*ar4[1] + epsilon[2]
ar4[3] = phi[1]*ar4[2] + phi[2]*ar4[1] + epsilon[3]
ar4[4] = phi[1]*ar4[3] + phi[2]*ar4[2] + phi[3]*ar4[1] + epsilon[4]

for t in 5:m
    ar4[t] = phi[1]*ar4[t-1] + phi[2]*ar4[t-2] + phi[3]*ar4[t-3] + phi[4]*ar4[t-4] + epsilon[t]
end


Y=primes_300[1:n]



M[1,1]=1
M[1,2]=-a_1
M[2,1]=-a_1
M[2,2]=a_1^2+1

M[1,3]=-a_2
M[3,1]=-a_2
M[2,3]=-a_1+a_1*a_2
M[3,2]=-a_1+a_1*a_2
M[3,3]=a_1^2+a_2^2+1

M[4,4]=a_1^2+a_2^2+a_3^2+1
M[4,3]=-a_1+a_1*a_2+a_3*a_2
M[4,2]=-a_2+a_1*a_3
M[4,1]=-a_3
M[3,4]=-a_1+a_1*a_2+a_3*a_2
M[2,4]=-a_2+a_1*a_3
M[1,4]=-a_3


M[n,n]=1
M[n,n-1]=-a_1
M[n-1,n]=-a_1
M[n-1,n-1]=a_1^2+1

M[n,n-2]=-a_2
M[n-2,n]=-a_2
M[n-1,n-2]=-a_1+a_1*a_2
M[n-2,n-1]=-a_1+a_1*a_2
M[n-2,n-2]=a_1^2+a_2^2+1

M[n-3,n-3]=a_1^2+a_2^2+a_3^2+1
M[n-3,n-2]=-a_1+a_1*a_2+a_3*a_2
M[n-3,n-1]=-a_2+a_1*a_3
M[n-3,n]=-a_3
M[n-2,n-3]=-a_1+a_1*a_2+a_3*a_2
M[n-1,n-3]=-a_2+a_1*a_3
M[n,n-3]=-a_3


D1[1,1]=0
D1[1,2]=-1
D1[2,1]=-1
D1[2,2]=2*a_1
D1[1,3]=0
D1[3,1]=0
D1[2,3]=-1+a_2
D1[3,2]=-1+a_2
D1[3,3]=2*a_1

D1[4,4]=2*a_1
D1[4,3]=-1+a_2
D1[4,2]=a_3
D1[4,1]=0
D1[3,4]=-1+a_2
D1[2,4]=a_3
D1[1,4]=0




D1[n,n]=0
D1[n,n-1]=-1
D1[n-1,n]=-1
D1[n-1,n-1]=2*a_1
D1[n,n-2]=0
D1[n-2,n]=0
D1[n-1,n-2]=-1+a_2
D1[n-2,n-1]=-1+a_2
D1[n-2,n-2]=2*a_1

D1[n-3,n-3]=2*a_1
D1[n-3,n-2]=-1+a_2
D1[n-3,n-1]=a_3
D1[n-3,n]=0
D1[n-2,n-3]=-1+a_2
D1[n-1,n-3]=a_3
D1[n,n-3]=0





D2[1,1]=0
D2[1,2]=0
D2[2,1]=0
D2[2,2]=0
D2[1,3]=-1
D2[3,1]=-1
D2[2,3]=a_1
D2[3,2]=a_1
D2[3,3]=2*a_2

D2[4,4]=2*a_2
D2[4,3]=a_1+a_3
D2[4,2]=-1
D2[4,1]=0
D2[3,4]=a_1+a_3
D2[2,4]=-1
D2[1,4]=0

D2[n,n]=0
D2[n,n-1]=0
D2[n-1,n]=0
D2[n-1,n-1]=0
D2[n,n-2]=-1
D2[n-2,n]=-1
D2[n-1,n-2]=a_1
D2[n-2,n-1]=a_1
D2[n-2,n-2]=2*a_2

D2[n-3,n-3]=2*a_2
D2[n-3,n-2]=a_1+a_3
D2[n-3,n-1]=-1
D2[n-3,n]=0
D2[n-2,n-3]=a_1+a_3
D2[n-1,n-3]=-1
D2[n,n-3]=0


D3[1,1]=0
D3[1,2]=0
D3[2,1]=0
D3[2,2]=0
D3[1,3]=0
D3[3,1]=0
D3[2,3]=0
D3[3,2]=0
D3[3,3]=0

D3[4,4]= 2*a_3
D3[4,3]=a_2
D3[4,2]=a_1
D3[4,1]=-1
D3[3,4]=a_2
D3[2,4]=a_1
D3[1,4]=-1

D3[n,n]=0
D3[n,n-1]=0
D3[n-1,n]=0
D3[n-1,n-1]=0
D3[n,n-2]=0
D3[n-2,n]=0
D3[n-1,n-2]=0
D3[n-2,n-1]=0
D3[n-2,n-2]=0

D3[n-3,n-3]=2*a_3
D3[n-3,n-2]=a_2
D3[n-3,n-1]=a_1
D3[n-3,n]=-1
D3[n-2,n-3]=a_2
D3[n-1,n-3]=a_1
D3[n,n-3]=-1


D4[1,1]=0
D4[1,2]=0
D4[2,1]=0
D4[2,2]=0
D4[1,3]=0
D4[3,1]=0
D4[2,3]=0
D4[3,2]=0
D4[3,3]=0

D4[4,4]=0
D4[4,3]=0
D4[4,2]=0
D4[4,1]=0
D4[3,4]=0
D4[2,4]=0
D4[1,4]=0

D4[n,n]=0
D4[n,n-1]=0
D4[n-1,n]=0
D4[n-1,n-1]=0
D4[n,n-2]=0
D4[n-2,n]=0
D4[n-1,n-2]=0
D4[n-2,n-1]=0
D4[n-2,n-2]=0

D4[n-3,n-3]=0 
D4[n-3,n-2]=0
D4[n-3,n-1]=0
D4[n-3,n]=0
D4[n-2,n-3]=0
D4[n-1,n-3]=0
D4[n,n-3]=0


T=[
    -a_1    1-a_2      -a_3    -a_4     0
    -a_2  -a_1-a_3     1-a_4      0     0
    -a_3   -a_2-a_4    -a_1       1    0
    -a_4    -a_3       -a_2    -a_1     1
      1     -a_1       -a_2    -a_3   -a_4
   
]

InvT=inv(T)

Gr=[0 
 0 
 0
 0
 1]

g3=InvT*Gr

g=[]
push!(g,g3[1])
push!(g,g3[2])
push!(g,g3[3])
push!(g,g3[4])
push!(g,g3[5])

g[1]
gd=(-a_1*((-a_1 - a_3)*(-a_2 - a_2*a_4 - a_1*(a_3 + a_1*a_4)) + (-a_1 - a_3*a_4 + (a_3 + a_1*a_4)*(-a_2 - a_4))*(-1 + a_4)) + ((a_1 + a_3)*(a_1*(-1 + a_4^2) - (a_2 + a_2*a_4)*a_3) - (a_1*(a_1 + a_3*a_4) + (-a_2 - a_4)*(a_2 + a_2*a_4))*a_2 + ((a_2 + a_4)*(-1 + a_4^2) - a_3*(a_1 + a_3*a_4))*(1 - a_4))*a_4 + (-(-a_2 - a_2*a_4 - a_1*(a_3 + a_1*a_4))*a_2 + (1 - (a_4^2) - (a_3 + a_1*a_4)*a_3)*(-1 + a_4))*(-1 + a_2) - ((a_1 + a_3)*(1 - (a_4^2) - (a_3 + a_1*a_4)*a_3) - (-a_1 - a_3*a_4 + (a_3 + a_1*a_4)*(-a_2 - a_4))*a_2)*a_3)


gm=g*0
gm[1]=simplify(g[1]*gd,expand=true)
gm[2]=simplify(g[2]*gd,expand=true)
gm[3]=simplify(g[3]*gd,expand=true)
gm[4]=simplify(g[4]*gd,expand=true)
gm[5]=simplify(g[5]*gd,expand=true)


gm[1]
gm[2]
gm[3]
gm[4]
gm[5]

for i in 6:n
    print(i)
    gtemp=a_1*gm[i-1]+a_2*gm[i-2]+a_3*gm[i-3]+a_4*gm[i-4]
    gtemp=simplify(gtemp,expand=true)
    push!(gm,gtemp)
end


Minv=zeros(Num, n, n)
for i in 1:n
    for j in 1:(n+1-i)
        Minv[j,j+i-1]=gm[i]
        Minv[j+i-1,j]=gm[i]
    end
end



deter=gd
Minvdet=Minv

contsol=[]
certsol=[]
noninvertible=[]
solf=[]

deter

deterfactor=(a_1-a_2+a_3-a_4+1)*(a_1+a_2+a_3+a_4-1)*(a_1^2*a_4-a_1*a_3*a_4+a_2*a_4^2-a_4^3+a_1*a_3-a_3^2+2*a_2*a_4-a_4^2+a_2+a_4+1)

factor1=(a_1-a_2+a_3-a_4+1)
factor2=(a_1+a_2+a_3+a_4-1)
factor3=(a_1^2*a_4-a_1*a_3*a_4+a_2*a_4^2-a_4^3+a_1*a_3-a_3^2+2*a_2*a_4-a_4^2+a_2+a_4+1)

#simplify(deter+factor1*factor2*factor3,expand=true)

tol=10^(-5)

Y=rand(-40:40,n)

for solutry in 1:100
    
    Y=rand(-100:100,n)
    @variables x y
    @variables a_1 a_2 a_3 a_4
    
    l1= tr(D1*Minvdet)*transpose(Y)*M*Y-n*transpose(Y)*D1*Y*deter
    l2= tr(D2*Minvdet)*transpose(Y)*M*Y-n*transpose(Y)*D2*Y*deter
    l3= tr(D3*Minvdet)*transpose(Y)*M*Y-n*transpose(Y)*D3*Y*deter
    l4= tr(D4*Minvdet)*transpose(Y)*M*Y-n*transpose(Y)*D4*Y*deter

    f1=simplify(l1,expand=true)
    f2=simplify(l2,expand=true)
    f3=simplify(l3,expand=true)
    f4=simplify(l4,expand=true)

    fa1 = string(f1)
    fa2 = string(f2)
    fa3 = string(f3)
    fa4 = string(f4)

    @var a0 a1 a2 a3 a4

    fa1_str = replace(fa1,"a_1" => "a1","a_2" => "a2","a_3" => "a3","a_4" => "a4")
    fa2_str = replace(fa2,"a_1" => "a1","a_2" => "a2","a_3" => "a3","a_4" => "a4")
    fa3_str = replace(fa3,"a_1" => "a1","a_2" => "a2","a_3" => "a3","a_4" => "a4")
    fa4_str = replace(fa4,"a_1" => "a1","a_2" => "a2","a_3" => "a3","a_4" => "a4")

    fa1_expr = eval(Meta.parse(fa1_str))
    fa2_expr = eval(Meta.parse(fa2_str))
    fa3_expr = eval(Meta.parse(fa3_str))
    fa4_expr = eval(Meta.parse(fa4_str))

    F = System([fa1_expr,fa2_expr,fa3_expr,fa4_expr])

    resultar4 = solve(F)

    cert=HomotopyContinuation.certify(F,resultar4)

    contf=0
    contne=0

    for i in 1:length(nonsingular(resultar4))

        if is_certified(certificates(cert)[i])
            b4=solution(nonsingular(resultar4)[i])[4]
            b3=solution(nonsingular(resultar4)[i])[3]
            b2=solution(nonsingular(resultar4)[i])[2]
            b1=solution(nonsingular(resultar4)[i])[1]
#            println(i," det= ",abs(bdgm), " cond1= ",abs(b1-b2+b3+1), " cond2= ", abs(b1+b2+b3-1)," cond3= ",abs(b1*b3-b3^2+b2+1) )
            #factor1=(a_1-a_2+a_3-a_4+1)
            #factor2=(a_1+a_2+a_3+a_4-1)
            #factor3=(a_1^2*a_4-a_1*a_3*a_4+a_2*a_4^2-a_4^3+a_1*a_3-a_3^2+2*a_2*a_4-a_4^2+a_2+a_4+1)
            if abs(b1-b2+b3-b4+1)<tol 
                contne=contne+1
            elseif abs(b1+b2+b3+b4-1)<tol
                contne=contne+1
            elseif abs(b1^2*b4-b1*b3*b4+b2*b4^2-b4^3+b1*b3-b3^2+2*b2*b4-b4^2+b2+b4+1)<tol
                contne=contne+1
            else 
                contf=contf+1
            end
        end 
    end


    println("Total: ",nresults(resultar4,only_nonsingular=true)," Cert: ",ndistinct_certified(cert), " Non_invertible: ",contne, " Final solutions: ", contf)
    push!(contsol,nresults(resultar4,only_nonsingular=true))
    push!(certsol,ndistinct_certified(cert))
    push!(noninvertible,contne)
    push!(solf,contf)


end

fa1_expr
fa2_expr
fa3_expr
fa4_expr

