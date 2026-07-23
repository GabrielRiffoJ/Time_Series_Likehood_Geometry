using HomotopyContinuation
using Symbolics
using LinearAlgebra


@variables a_1 a_2 a_3 a_4 a_5 a_6 a_7 a_8 a_9 a_10

n=20

M = zeros(Num, n, n)
    
for i in 1:n
    M[i, i] =a_6^2+a_5^2+a_4^2+a_3^2+a_2^2+a_1^2+1
end

for i in 1:(n-1)
    M[i, i+1] = -a_1+a_1*a_2+a_2*a_3+a_3*a_4+a_4*a_5+a_5*a_6
    M[i+1, i] = -a_1+a_1*a_2+a_2*a_3+a_3*a_4+a_4*a_5+a_5*a_6
end

for i in 1:(n-2)
    M[i, i+2] = -a_2+a_1*a_3+a_2*a_4+a_3*a_5+a_4*a_6
    M[i+2, i] = -a_2+a_1*a_3+a_2*a_4+a_3*a_5+a_4*a_6
end

for i in 1:(n-3)
    M[i, i+3] = -a_3+a_1*a_4+a_2*a_5+a_3*a_6
    M[i+3, i] = -a_3+a_1*a_4+a_2*a_5+a_3*a_6
end

for i in 1:(n-4)
    M[i, i+4] = -a_4+a_1*a_5+a_4*a_6
    M[i+4, i] = -a_4+a_1*a_5+a_4*a_6
end

for i in 1:(n-5)
    M[i, i+5] = -a_5+a_5*a_6
    M[i+5, i] = -a_5+a_5*a_6
end

for i in 1:(n-6)
    M[i, i+6] = -a_6
    M[i+6, i] = -a_6
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

for i in 1:(n-4)
    D1[i, i+4] = a_5
    D1[i+4, i] = a_5
end

for i in 1:(n-5)
    D1[i, i+5] = a_6
    D1[i+5, i] = a_6
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

for i in 1:(n-3)
    D2[i, i+3] = a_5
    D2[i+3, i] = a_5
end

for i in 1:(n-4)
    D2[i, i+4] = a_6
    D2[i+4, i] = a_6
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
    D3[i, i+2] = a_1+a_5
    D3[i+2, i] = a_1+a_5
end

for i in 1:(n-3)
    D3[i, i+3] = -1+a_6
    D3[i+3, i] = -1+a_6
end

D4 = zeros(Num, n, n)

for i in 1:(n)
    D4[i, i] = 2*a_4
end


for i in 1:(n-1)
    D4[i, i+1] = a_3+a_5
    D4[i+1, i] = a_3+a_5
end

for i in 1:(n-2)
    D4[i, i+2] = a_2+a_6
    D4[i+2, i] = a_2+a_6
end

for i in 1:(n-3)
    D4[i, i+3] = a_1
    D4[i+3, i] = a_1
end

for i in 1:(n-4)
    D4[i, i+4] = -1
    D4[i+4, i] = -1
end


D5 = zeros(Num, n, n)

for i in 1:(n)
    D5[i, i] = 2*a_5
end


for i in 1:(n-1)
    D5[i, i+1] = a_4+a_6
    D5[i+1, i] = a_4+a_6
end

for i in 1:(n-2)
    D5[i, i+2] = a_3
    D5[i+2, i] = a_3
end

for i in 1:(n-3)
    D5[i, i+3] = a_2
    D5[i+3, i] = a_2
end

for i in 1:(n-4)
    D5[i, i+4] = a_1
    D5[i+4, i] = a_1
end

for i in 1:(n-5)
    D5[i, i+5] = -1
    D5[i+5, i] = -1
end

D6 = zeros(Num, n, n)

for i in 1:(n)
    D6[i, i] = 2*a_6
end


for i in 1:(n-1)
    D6[i, i+1] = a_5
    D6[i+1, i] = a_5
end

for i in 1:(n-2)
    D6[i, i+2] = a_4
    D6[i+2, i] = a_4 
end

for i in 1:(n-3)
    D6[i, i+3] = a_3
    D6[i+3, i] = a_3
end

for i in 1:(n-4)
    D6[i, i+4] = a_2
    D6[i+4, i] = a_2
end

for i in 1:(n-5)
    D6[i, i+5] = a_1
    D6[i+5, i] = a_1
end

for i in 1:(n-6)
    D6[i, i+6] = -1
    D6[i+6, i] = -1
end



# === MATRIX M (AR6) ===

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

M[5,5]=a_1^2+a_2^2+a_3^2+a_4^2+1
M[5,4]=-a_1+a_1*a_2+a_3*a_2+a_4*a_3
M[5,3]=-a_2+a_1*a_3+a_2*a_4
M[5,2]=-a_3+a_1*a_4
M[5,1]=-a_4
M[4,5]=-a_1+a_1*a_2+a_3*a_2+a_4*a_3
M[3,5]=-a_2+a_1*a_3+a_2*a_4
M[2,5]=-a_3+a_1*a_4
M[1,5]=-a_4

M[6,6]=a_1^2+a_2^2+a_3^2+a_4^2+a_5^2+1
M[6,5]=-a_1+a_1*a_2+a_3*a_2+a_4*a_3+a_5*a_4
M[6,4]=-a_2+a_1*a_3+a_2*a_4+a_3*a_5
M[6,3]=-a_3+a_1*a_4+a_2*a_5
M[6,2]=-a_4+a_1*a_5
M[6,1]=-a_5
M[5,6]=-a_1+a_1*a_2+a_3*a_2+a_4*a_3+a_5*a_4
M[4,6]=-a_2+a_1*a_3+a_2*a_4+a_3*a_5
M[3,6]=-a_3+a_1*a_4+a_2*a_5
M[2,6]=-a_4+a_1*a_5
M[1,6]=-a_5


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

M[n-4,n-4]=a_1^2+a_2^2+a_3^2+a_4^2+1
M[n-4,n-3]=-a_1+a_1*a_2+a_3*a_2+a_4*a_3
M[n-4,n-2]=-a_2+a_1*a_3+a_2*a_4
M[n-4,n-1]=-a_3+a_1*a_4
M[n-4,n]=-a_4
M[n-3,n-4]=-a_1+a_1*a_2+a_3*a_2+a_4*a_3
M[n-2,n-4]=-a_2+a_1*a_3+a_2*a_4
M[n-1,n-4]=-a_3+a_1*a_4
M[n,n-4]=-a_4

M[n-5,n-5]=a_1^2+a_2^2+a_3^2+a_4^2+a_5^2+1
M[n-5,n-4]=-a_1+a_1*a_2+a_3*a_2+a_4*a_3+a_5*a_4
M[n-5,n-3]=-a_2+a_1*a_3+a_2*a_4+a_3*a_5
M[n-5,n-2]=-a_3+a_1*a_4+a_2*a_5
M[n-5,n-1]=-a_4+a_1*a_5
M[n-5,n]=-a_5
M[n-4,n-5]=-a_1+a_1*a_2+a_3*a_2+a_4*a_3+a_5*a_4
M[n-3,n-5]=-a_2+a_1*a_3+a_2*a_4+a_3*a_5
M[n-2,n-5]=-a_3+a_1*a_4+a_2*a_5
M[n-1,n-5]=-a_4+a_1*a_5
M[n,n-5]=-a_5


# === D1 ===

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

D1[5,5]=2*a_1
D1[5,4]=-1+a_2
D1[5,3]=a_3
D1[5,2]=a_4
D1[5,1]=0
D1[4,5]=-1+a_2
D1[3,5]=a_3
D1[2,5]=a_4
D1[1,5]=0

D1[6,6]=2*a_1
D1[6,5]=-1+a_2
D1[6,4]=a_3
D1[6,3]=a_4
D1[6,2]=a_5
D1[6,1]=0
D1[5,6]=-1+a_2
D1[4,6]=a_3
D1[3,6]=a_4
D1[2,6]=a_5
D1[1,6]=0

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

D1[n-4,n-4]=2*a_1
D1[n-4,n-3]=-1+a_2
D1[n-4,n-2]=a_3
D1[n-4,n-1]=a_4
D1[n-4,n]=0
D1[n-3,n-4]=-1+a_2
D1[n-2,n-4]=a_3
D1[n-1,n-4]=a_4
D1[n,n-4]=0

D1[n-5,n-5]=2*a_1
D1[n-5,n-4]=-1+a_2
D1[n-5,n-3]=a_3
D1[n-5,n-2]=a_4
D1[n-5,n-1]=a_5
D1[n-5,n]=0
D1[n-4,n-5]=-1+a_2
D1[n-3,n-5]=a_3
D1[n-2,n-5]=a_4
D1[n-1,n-5]=a_5
D1[n,n-5]=0


# === D2 ===

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

D2[5,5]=2*a_2
D2[5,4]=a_1+a_3
D2[5,3]=-1+a_4
D2[5,2]=0
D2[5,1]=0
D2[4,5]=a_1+a_3
D2[3,5]=-1+a_4
D2[2,5]=0
D2[1,5]=0

D2[6,6]=2*a_2
D2[6,5]=a_1+a_3
D2[6,4]=-1+a_4
D2[6,3]=a_5
D2[6,2]=0
D2[6,1]=0
D2[5,6]=a_1+a_3
D2[4,6]=-1+a_4
D2[3,6]=a_5
D2[2,6]=0
D2[1,6]=0

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

D2[n-4,n-4]=2*a_2
D2[n-4,n-3]=a_1+a_3
D2[n-4,n-2]=-1+a_4
D2[n-4,n-1]=0
D2[n-4,n]=0
D2[n-3,n-4]=a_1+a_3
D2[n-2,n-4]=-1+a_4
D2[n-1,n-4]=0
D2[n,n-4]=0

D2[n-5,n-5]=2*a_2
D2[n-5,n-4]=a_1+a_3
D2[n-5,n-3]=-1+a_4
D2[n-5,n-2]=a_5
D2[n-5,n-1]=0
D2[n-5,n]=0
D2[n-4,n-5]=a_1+a_3
D2[n-3,n-5]=-1+a_4
D2[n-2,n-5]=a_5
D2[n-1,n-5]=0
D2[n,n-5]=0


# === D3 ===


D3[1,1]=0
D3[1,2]=0
D3[2,1]=0
D3[2,2]=0
D3[1,3]=0
D3[3,1]=0
D3[2,3]=0
D3[3,2]=0
D3[3,3]=0

D3[4,4]=2*a_3
D3[4,3]=a_2
D3[4,2]=a_1
D3[4,1]=-1
D3[3,4]=a_2
D3[2,4]=a_1
D3[1,4]=-1

D3[5,5]=2*a_3
D3[5,4]=a_2+a_4
D3[5,3]=a_1
D3[5,2]=-1
D3[5,1]=0
D3[4,5]=a_2+a_4
D3[3,5]=a_1
D3[2,5]=-1
D3[1,5]=0

D3[6,6]=2*a_3
D3[6,5]=a_2+a_4
D3[6,4]=a_1+a_5
D3[6,3]=-1
D3[6,2]=0
D3[6,1]=0
D3[5,6]=a_2+a_4
D3[4,6]=a_1+a_5
D3[3,6]=-1
D3[2,6]=0
D3[1,6]=0


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

D3[n-4,n-4]=2*a_3
D3[n-4,n-3]=a_2+a_4
D3[n-4,n-2]=a_1
D3[n-4,n-1]=-1
D3[n-4,n]=0
D3[n-3,n-4]=a_2+a_4
D3[n-2,n-4]=a_1
D3[n-1,n-4]=-1
D3[n,n-4]=0

D3[n-5,n-5]=2*a_3
D3[n-5,n-4]=a_2+a_4
D3[n-5,n-3]=a_1+a_5
D3[n-5,n-2]=-1
D3[n-5,n-1]=0
D3[n-5,n]=0
D3[n-4,n-5]=a_2+a_4
D3[n-3,n-5]=a_1+a_5
D3[n-2,n-5]=-1
D3[n-1,n-5]=0
D3[n,n-5]=0


# === D4 ===


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

D4[5,5]=2*a_4
D4[5,4]=a_3
D4[5,3]=a_2
D4[5,2]=a_1
D4[5,1]=-1
D4[4,5]=a_3
D4[3,5]=a_2
D4[2,5]=a_1
D4[1,5]=-1

D4[6,6]=2*a_4
D4[6,5]=a_3+a_5
D4[6,4]=a_2
D4[6,3]=a_1
D4[6,2]=-1
D4[6,1]=0
D4[5,6]=a_3+a_5
D4[4,6]=a_2
D4[3,6]=a_1
D4[2,6]=-1
D4[1,6]=0


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

D4[n-4,n-4]=2*a_4
D4[n-4,n-3]=a_3
D4[n-4,n-2]=a_2
D4[n-4,n-1]=a_1
D4[n-4,n]=-1
D4[n-3,n-4]=a_3
D4[n-2,n-4]=a_2
D4[n-1,n-4]=a_1
D4[n,n-4]=-1

D4[n-5,n-5]=2*a_4
D4[n-5,n-4]=a_3+a_5
D4[n-5,n-3]=a_2
D4[n-5,n-2]=a_1
D4[n-5,n-1]=-1
D4[n-5,n]=0
D4[n-4,n-5]=a_3+a_5
D4[n-3,n-5]=a_2
D4[n-2,n-5]=a_1
D4[n-1,n-5]=-1
D4[n,n-5]=0


# === D5 ===

D5[1,1]=0
D5[1,2]=0
D5[2,1]=0
D5[2,2]=0
D5[1,3]=0
D5[3,1]=0
D5[2,3]=0
D5[3,2]=0
D5[3,3]=0

D5[4,4]=0
D5[4,3]=0
D5[4,2]=0
D5[4,1]=0
D5[3,4]=0
D5[2,4]=0
D5[1,4]=0

D5[5,5]=0
D5[5,4]=0
D5[5,3]=0
D5[5,2]=0
D5[5,1]=0
D5[4,5]=0
D5[3,5]=0
D5[2,5]=0
D5[1,5]=0

D5[6,6]=2*a_5
D5[6,5]=a_4
D5[6,4]=a_3
D5[6,3]=a_2
D5[6,2]=a_1
D5[6,1]=-1
D5[5,6]=a_4
D5[4,6]=a_3
D5[3,6]=a_2
D5[2,6]=a_1
D5[1,6]=-1

D5[n,n]=0
D5[n,n-1]=0
D5[n-1,n]=0
D5[n-1,n-1]=0
D5[n,n-2]=0
D5[n-2,n]=0
D5[n-1,n-2]=0
D5[n-2,n-1]=0
D5[n-2,n-2]=0

D5[n-3,n-3]=0
D5[n-3,n-2]=0
D5[n-3,n-1]=0
D5[n-3,n]=0
D5[n-2,n-3]=0
D5[n-1,n-3]=0
D5[n,n-3]=0

D5[n-4,n-4]=0
D5[n-4,n-3]=0
D5[n-4,n-2]=0
D5[n-4,n-1]=0
D5[n-4,n]=0
D5[n-3,n-4]=0
D5[n-2,n-4]=0
D5[n-1,n-4]=0
D5[n,n-4]=0

D5[n-5,n-5]=2*a_5
D5[n-5,n-4]=a_4
D5[n-5,n-3]=a_3
D5[n-5,n-2]=a_2
D5[n-5,n-1]=a_1
D5[n-5,n]=-1
D5[n-4,n-5]=a_4
D5[n-3,n-5]=a_3
D5[n-2,n-5]=a_2
D5[n-1,n-5]=a_1
D5[n,n-5]=-1


# === D6 ===


D6[1,1]=0
D6[1,2]=0
D6[2,1]=0
D6[2,2]=0
D6[1,3]=0
D6[3,1]=0
D6[2,3]=0
D6[3,2]=0
D6[3,3]=0

D6[4,4]=0
D6[4,3]=0
D6[4,2]=0
D6[4,1]=0
D6[3,4]=0
D6[2,4]=0
D6[1,4]=0

D6[5,5]=0
D6[5,4]=0
D6[5,3]=0
D6[5,2]=0
D6[5,1]=0
D6[4,5]=0
D6[3,5]=0
D6[2,5]=0
D6[1,5]=0

D6[6,6]=0
D6[6,5]=0
D6[6,4]=0
D6[6,3]=0
D6[6,2]=0
D6[6,1]=0
D6[5,6]=0
D6[4,6]=0
D6[3,6]=0
D6[2,6]=0
D6[1,6]=0


D6[n,n]=0
D6[n,n-1]=0
D6[n-1,n]=0
D6[n-1,n-1]=0
D6[n,n-2]=0
D6[n-2,n]=0
D6[n-1,n-2]=0
D6[n-2,n-1]=0
D6[n-2,n-2]=0

D6[n-3,n-3]=0
D6[n-3,n-2]=0
D6[n-3,n-1]=0
D6[n-3,n]=0
D6[n-2,n-3]=0
D6[n-1,n-3]=0
D6[n,n-3]=0

D6[n-4,n-4]=0
D6[n-4,n-3]=0
D6[n-4,n-2]=0
D6[n-4,n-1]=0
D6[n-4,n]=0
D6[n-3,n-4]=0
D6[n-2,n-4]=0
D6[n-1,n-4]=0
D6[n,n-4]=0

D6[n-5,n-5]=0
D6[n-5,n-4]=0
D6[n-5,n-3]=0
D6[n-5,n-2]=0
D6[n-5,n-1]=0
D6[n-5,n]=0
D6[n-4,n-5]=0
D6[n-3,n-5]=0
D6[n-2,n-5]=0
D6[n-1,n-5]=0
D6[n,n-5]=0


D6[n-6,n-6]=0
D6[n-6,n-5]=0
D6[n-6,n-4]=0
D6[n-6,n-3]=0
D6[n-6,n-2]=0
D6[n-6,n-1]=0
D6[n-6,n]=0
D6[n-5,n-6]=0
D6[n-4,n-6]=0
D6[n-3,n-6]=0
D6[n-2,n-6]=0
D6[n-1,n-6]=0
D6[n,n-6]=0

using Random
using Primes


prim = Primes.primes(1, 2000)      
primes_300 = prim[1:300]             

inv_primes = 1.0 ./ primes_300

m = 300
phi = [0.4, 0.3, 0.2, 0.1, 0.05, 0.01]
sigma = 1.0
Random.seed!(1234)

ar6 = zeros(m)
epsilon = sigma .* randn(m)

ar6[1] = epsilon[1]
ar6[2] = phi[1]*ar6[1] + epsilon[2]
ar6[3] = phi[1]*ar6[2] + phi[2]*ar6[1] + epsilon[3]
ar6[4] = phi[1]*ar6[3] + phi[2]*ar6[2] + phi[3]*ar6[1] + epsilon[4]
ar6[5] = phi[1]*ar6[4] + phi[2]*ar6[3] + phi[3]*ar6[2] + phi[4]*ar6[1] + epsilon[5]
ar6[6] = phi[1]*ar6[5] + phi[2]*ar6[4] + phi[3]*ar6[3] + phi[4]*ar6[2] + phi[5]*ar6[1] + epsilon[6]

for t in 7:m
    ar6[t] = phi[1]*ar6[t-1] + phi[2]*ar6[t-2] + phi[3]*ar6[t-3] + phi[4]*ar6[t-4] + phi[5]*ar6[t-5] + phi[6]*ar6[t-6] + epsilon[t]
end


Y=primes_300[1:n]


deter=gd
Minvdet=Minv


#for i in 1:n
#    for j in 1:n
#        Minvdet[i,j]=simplify(Minv[i,j]*deter,expand=true)
#    end
#end

@variables x y
@variables a_1 a_2

trace1=0
trace2=0
trace3=0
trace4=0
trace5=0
trace6=0

for i in 1:n
    for j in 1:n
        trace1=trace1+D1[i,j]*Minvdet[j,i]
        trace2=trace2+D2[i,j]*Minvdet[j,i]
        trace3=trace3+D3[i,j]*Minvdet[j,i]
        trace4=trace4+D4[i,j]*Minvdet[j,i]
        trace5=trace5+D5[i,j]*Minvdet[j,i]
        trace6=trace6+D6[i,j]*Minvdet[j,i]
    end
end

trace1_exp=simplify(trace1,expand=true)
trace2_exp=simplify(trace2,expand=true)
trace3_exp=simplify(trace3,expand=true)
trace4_exp=simplify(trace4,expand=true)
trace5_exp=simplify(trace5,expand=true)
trace6_exp=simplify(trace6,expand=true)


contsol=[]
certsol=[]
noninvertible=[]
solf=[]

for solutry in 1:100
    
    Y=rand(-40:40,n)

    l1= trace1_exp*transpose(Y)*M*Y-n*transpose(Y)*D1*Y*deter
    l2= trace2_exp*transpose(Y)*M*Y-n*transpose(Y)*D2*Y*deter
    l3= trace3_exp*transpose(Y)*M*Y-n*transpose(Y)*D3*Y*deter
    l4= trace4_exp*transpose(Y)*M*Y-n*transpose(Y)*D4*Y*deter
    l5= trace5_exp*transpose(Y)*M*Y-n*transpose(Y)*D5*Y*deter
    l6= trace6_exp*transpose(Y)*M*Y-n*transpose(Y)*D6*Y*deter

    f1=simplify(l1,expand=true)
    f2=simplify(l2,expand=true)
    f3=simplify(l3,expand=true)
    f4=simplify(l4,expand=true)
    f5=simplify(l5,expand=true)
    f6=simplify(l6,expand=true)

    fa1 = string(f1)
    fa2 = string(f2)
    fa3 = string(f3)
    fa4 = string(f4)
    fa5 = string(f5)
    fa6 = string(f6)

    @var a0 a1 a2 a3 a4 a5 a6

    fa1_str = replace(fa1,"a_1" => "a1","a_2" => "a2","a_3" => "a3","a_4" => "a4","a_5" => "a5","a_6" => "a6")
    fa2_str = replace(fa2,"a_1" => "a1","a_2" => "a2","a_3" => "a3","a_4" => "a4","a_5" => "a5","a_6" => "a6")
    fa3_str = replace(fa3,"a_1" => "a1","a_2" => "a2","a_3" => "a3","a_4" => "a4","a_5" => "a5","a_6" => "a6")
    fa4_str = replace(fa4,"a_1" => "a1","a_2" => "a2","a_3" => "a3","a_4" => "a4","a_5" => "a5","a_6" => "a6")
    fa5_str = replace(fa5,"a_1" => "a1","a_2" => "a2","a_3" => "a3","a_4" => "a4","a_5" => "a5","a_6" => "a6")
    fa6_str = replace(fa6,"a_1" => "a1","a_2" => "a2","a_3" => "a3","a_4" => "a4","a_5" => "a5","a_6" => "a6")

    fa1_expr = eval(Meta.parse(fa1_str))
    fa2_expr = eval(Meta.parse(fa2_str))
    fa3_expr = eval(Meta.parse(fa3_str))
    fa4_expr = eval(Meta.parse(fa4_str))
    fa5_expr = eval(Meta.parse(fa5_str))
    fa6_expr = eval(Meta.parse(fa6_str))

    F = System([fa1_expr,fa2_expr,fa3_expr,fa4_expr,fa5_expr,fa6_expr])
    
    resultar6 = HomotopyContinuation.solve(F)


    cert=HomotopyContinuation.certify(F,resultar6)
    
    i=1
    
    contf=0
    contne=0

    for i in 1:length(nonsingular(resultar6))

        if is_certified(certificates(cert)[i])
                
                b6=solution(nonsingular(resultar6)[i])[6]
                b5=solution(nonsingular(resultar6)[i])[5]
                b4=solution(nonsingular(resultar6)[i])[4]
                b3=solution(nonsingular(resultar6)[i])[3]
                b2=solution(nonsingular(resultar6)[i])[2]
                b1=solution(nonsingular(resultar6)[i])[1]
                expr_sub=abs(substitute(deter, Dict(a_1 => b1,a_2 => b2,a_3 => b3, a_4 => b4,a_5 => b5,a_6 => b6)))
                
                if expr_sub <1e-8
                    contne=contne+1
                else 
                    contf=contf+1
                end

        end

    end
    println("Total: ",nresults(resultar6,only_nonsingular=true)," Cert: ",ndistinct_certified(cert), " Non_stationary: ",contne, " Final solutions: ", contf)
    push!(contsol,nresults(resultar6,only_nonsingular=true))
    push!(certsol,ndistinct_certified(cert))
    push!(noninvertible,contne)
    push!(solf,contf)
end

