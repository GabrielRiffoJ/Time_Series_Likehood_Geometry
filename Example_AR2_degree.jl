using HomotopyContinuation
using Symbolics
using LinearAlgebra



@variables a_1 a_2 a_3 a_4 a_5 a_6 a_7 a_8 a_9 a_10

n=10

M = zeros(Num, n, n)
    

for i in 1:n
    M[i, i] = a_2^2+a_1^2+1
end

for i in 1:(n-1)
    M[i, i+1] = -a_1+a_1*a_2
    M[i+1, i] = -a_1+a_1*a_2
end

for i in 1:(n-2)
    M[i, i+2] = -a_2
    M[i+2, i] = -a_2
end


D1 = zeros(Num, n, n)

for i in 1:(n)
    D1[i, i] = 2*a_1
end


for i in 1:(n-1)
    D1[i, i+1] = -1+a_2
    D1[i+1, i] = -1+a_2
end

D2 = zeros(Num, n, n)

for i in 1:(n)
    D2[i, i] = 2*a_2
end


for i in 1:(n-1)
    D2[i, i+1] = a_1
    D2[i+1, i] = a_1
end

for i in 1:(n-2)
    D2[i, i+2] = -1
    D2[i+2, i] = -1
end


YI=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73]
YF=[1/2,1/3,1/5,1/7,1/11,1/13,1/17,1/19,1/23,1/29,1/31,1/37]
Y=YI[2:(n+1)]

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


D2[1,1]=0
D2[1,2]=0
D2[2,1]=0
D2[2,2]=0

D2[n,n]=0
D2[n-1,n]=0
D2[n,n-1]=0
D2[n-1,n-1]=0



#M=[
#1	   -a_1	                  -a_2	              0	                    0	        0	        0
#-a_1	a_1^2+1	         -a_1+a_1*a_2	        -a_2	                0	        0	        0
#-a_2    -a_1+a_1*a_2	    a_2^2+a_1^2+1	    -a_1+a_1*a_2	      -a_2	        0	        0
# 0	   -a_2	                 -a_1+a_1*a_2	    a_1^2+1+a_2^2	-a_1+a_1*a_2	   -a_2	        0
# 0	    0	                 -a_2	            -a_1+a_1*a_2	a_1^2+1+a_2^2	-a_1+a_1*a_2  -a_2
# 0	    0	                    0	             -a_2	          -a_1+a_1*a_2	  a_1^2+1	  -a_1
# 0	    0	                    0	              0	                  -a_2	      -a_1	        1
#]

Minv=inv(M)
deter=det(M)
Minvdet=Minv*deter

for i in 1:n
    for j in 1:n
        Minvdet[i,j]=simplify(Minv[i,j]*deter,expand=true)
    end
end

#Sinvdet[1,5]
#simplify(Sinv[1,1],expand=true)


#D1=[
#0	   -1	                    0	              0	                0	        0	        0
#-1	   2*a_1	              -1+a_2	          0	                0	        0	        0
#0     -1+a_2	              2*a_1	            -1+a_2	            0	        0	        0
# 0	   0	                 -1+a_2	            2*a_1	         -1+a_2	        0           0
# 0	    0	                    0 	            -1+a_2	          2*a_1^2	 -1+a_2         0
# 0	    0	                    0	              0	              -1+a_2	  2*a_1	       -1
# 0	    0	                    0	              0	                0	      -1	        0 
#]


#D2=[
#0	    0	                  -1	              0	                    0	        0	        0
#0	    0	                 -a_1	             -1	                    0	        0	        0
#-1     a_1	                 2*a_2	             a_1	               -1	        0	        0
# 0	   -1	                  a_1	            2*a_2       	       a_1         -1	        0
# 0	    0	                   -1	             a_1	              2*a_2	       a_1          -1
# 0	    0	                    0	             -1	                   a_1	        0	        0
# 0	    0	                    0	              0	                   -1	        0	        0
#]

using Random
Random.seed!(1234)


@variables x y
@variables a_1 a_2

Y=YF[1:(n)]
    
l1= tr(D1*Minvdet)*transpose(Y)*M*Y-n*transpose(Y)*D1*Y*deter
l2= tr(D2*Minvdet)*transpose(Y)*M*Y-n*transpose(Y)*D2*Y*deter

f1=simplify(l1,expand=true)
f2=simplify(l2,expand=true)

fa1 = string(f1)
fa2 = string(f2)

@var a0 a1 a2

fa1_str = replace(fa1,"a_1" => "a1","a_2" => "a2")
fa2_str = replace(fa2,"a_1" => "a1","a_2" => "a2")

fa1_expr = eval(Meta.parse(fa1_str))
fa2_expr = eval(Meta.parse(fa2_str))

F = System([fa1_expr,fa2_expr])

resultcl2 = HomotopyContinuation.solve(F)

