
R = QQ[a_0, a_1, a_2, a_3,y_1, y_2, y_3, y_4]

--R = (frac QQ[y_1, y_2, y_3, y_4])[a_0, a_1, a_2, a_3]

K = frac(R)


M = matrix{{a_0, a_1, a_2, a_3},
           {a_1, a_0, a_1, a_2},
           {a_2, a_1, a_0, a_1},
           {a_3, a_2, a_1, a_0}}


M0 = matrix{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}
M1 = matrix{{0,1,0,0},{1,0,1,0},{0,1,0,1},{0,0,1,0}}
M2 = matrix{{0,0,1,0},{0,0,0,1},{1,0,0,0},{0,1,0,0}}
M3 = matrix{{0,0,0,1},{0,0,0,0},{0,0,0,0},{1,0,0,0}}

deter = det(M)


Minv = inverse(M)


Minvdet = Minv * deter


Y = matrix{{y_1},{y_2},{y_3},{y_1+y_2-y_3}}

--- 2,3,5,0
--Degenerate case 0
--Y = matrix{{-2*y_2-y_3},{y_2},{y_3},{-y_2-2*y_3}}

--Y = matrix{{-y_2},{y_2},{-y_2},{y_2}}


l0 = trace(Minvdet * M0) * deter - transpose(Y) * Minvdet * M0 * Minvdet * Y
l1 = trace(Minvdet * M1) * deter - transpose(Y) * Minvdet * M1 * Minvdet * Y
l2 = trace(Minvdet * M2) * deter - transpose(Y) * Minvdet * M2 * Minvdet * Y
l3 = trace(Minvdet * M3) * deter - transpose(Y) * Minvdet * M3 * Minvdet * Y


f0 = l0_0_0
f1 = l1_0_0
f2 = l2_0_0
f3 = l3_0_0

f0R = numerator(f0)
f1R = numerator(f1)
f2R = numerator(f2)
f3R = numerator(f3)

ring f0R
ring f1R
ring f2R
ring f3R

---E1 = 2*y_1^3+3*y_1^2*y_2-3*y_1^2*y_3+18*y_1*y_2*y_3-9*y_2^3+21*y_2^2*y_3-8*y_3^3
--FE1=numerator(E1)
--I=ideal(f0R,f1R,f2R,f3R,FE1)

I=ideal(f0R,f1R,f2R,f3R)

detR=numerator(deter)


Mn=minimalPrimes (ideal(detR))

fac1=Mn_0_0
fac2=Mn_1_0



J1 = saturate(I, fac1)


J2 = saturate(J1, fac2)



Jtrim=trim(J2)


--RJ=radical(J2)

RJ2=RJ


--MPRJ= minimalPrimes RJ


S= QQ[y_1, y_2, y_3,y_4,a_3,a_2, a_1, a_0, MonomialOrder => Lex]

Say= QQ[y_1, y_2, y_3,y_4,a_0,a_1, a_2, a_3, MonomialOrder => Lex]

Say= QQ[a_2,a_3,y_3,y_4,MonomialOrder => Lex]

S3= QQ[y_3,y_4,a_3][a_2]
S1= QQ[y_3,y_4,a_3,a_2][a_1]
S0= QQ[y_3,y_4,a_3,a_2,a_1][a_0]

Sgb=QQ[y_1,y_2,y_3,y_4][a_3,a_2,a_1,a_0, MonomialOrder => Lex]

RJ2= substitute(RJ, S)

GRJ2 = gens gb RJ2

RJ2




Eja3=(729*y_3^6-891*y_3^4*y_4^2+171*y_3^2*y_4^4-9*y_4^6-972*y_3^5+1458*y_3^4*y_4-1512*y_3^3*y_
       4^2+2268*y_3^2*y_4^3+180*y_3*y_4^4-270*y_4^5-22599*y_3^4+19440*y_3^3*y_4-342*y_3^2*y_4^2+1296*y
       _3*y_4^3-2403*y_4^4+18360*y_3^3-11988*y_3^2*y_4-22824*y_3*y_4^2-756*y_4^3+175779*y_3^2-180144*y
       _3*y_4+74817*y_4^2-130284*y_3+195426*y_4-266805)

Eja3c2=(576*y_3^7*y_4-480*y_3^5*y_4^3+48*y_3^3*
       y_4^5-5022*y_3^7-3054*y_3^6*y_4+5994*y_3^5*y_4^2+1986*y_3^4*y_4^3+522*y_3^3*y_4^4-498*y_3^2*y_4
       ^5-54*y_3*y_4^6+30*y_4^7+8964*y_3^6-22902*y_3^5*y_4+14688*y_3^4*y_4^2-11892*y_3^3*y_4^3-6012*y_
       3^2*y_4^4-366*y_3*y_4^5+792*y_4^6+153522*y_3^5-33918*y_3^4*y_4-60840*y_3^3*y_4^2-18468*y_3^2*y_
       4^3+13410*y_3*y_4^4+5850*y_4^5-210600*y_3^4+213852*y_3^3*y_4+67032*y_3^2*y_4^2+135012*y_3*y_4^3
       -8928*y_4^4-1192698*y_3^3+471126*y_3^2*y_4+68670*y_3*y_4^2-230382*y_4^3+1654884*y_3^2-1746438*y
       _3*y_4-289944*y_4^2+1647270*y_3+1652694*y_4-1627776)

Eja3c1=(125*y_3^10+53*y_3^8*y_4^2-85*y_3^6*
       y_4^4+3*y_3^4*y_4^6-900*y_3^9-822*y_3^8*y_4-732*y_3^7*y_4^2+822*y_3^6*y_4^3+764*y_3^5*y_4^4+306
       *y_3^4*y_4^5-92*y_3^3*y_4^6-18*y_3^2*y_4^7+4767*y_3^8+8940*y_3^7*y_4-1644*y_3^6*y_4^2-6672*y_3^
       5*y_4^3-9172*y_3^4*y_4^4-348*y_3^3*y_4^5+742*y_3^2*y_4^6-y_4^8+19480*y_3^7-7284*y_3^6*y_4-12868
       *y_3^5*y_4^2+39648*y_3^4*y_4^3+9788*y_3^3*y_4^4+13638*y_3^2*y_4^5-1048*y_3*y_4^6-534*y_4^7-
       235722*y_3^6-1236*y_3^5*y_4-81556*y_3^4*y_4^2+202464*y_3^3*y_4^3-29515*y_3^2*y_4^4-17076*y_3*y_
       4^5-8169*y_4^6+3360*y_3^5+631632*y_3^4*y_4-600896*y_3^3*y_4^2-450192*y_3^2*y_4^3-123844*y_3*y_4
       ^4+3570*y_4^5+1865066*y_3^4-749964*y_3^3*y_4+1266732*y_3^2*y_4^2-41424*y_3*y_4^3+397860*y_4^4-
       2133480*y_3^3-1938564*y_3^2*y_4+3222728*y_3*y_4^2+96486*y_4^3+541581*y_3^2-4305420*y_3*y_4-
       4974209*y_4^2+1739988*y_3+10653390*y_4-7796169)



Eja3c0=      (8*y_3^9*y_4^3-8*y_3^7*y_4^5-250*y_3^10*y_4-
       198*y_3^9*y_4^2+66*y_3^8*y_4^3+72*y_3^7*y_4^4+58*y_3^6*y_4^5+30*y_3^5*y_4^6-2*y_3^4*y_4^7+2696*
       y_3^9*y_4+276*y_3^8*y_4^2-568*y_3^7*y_4^3-180*y_3^6*y_4^4-1076*y_3^5*y_4^5-300*y_3^4*y_4^6+76*y
       _3^3*y_4^7+12*y_3^2*y_4^8+936*y_3^9-6286*y_3^8*y_4-3804*y_3^7*y_4^2+1352*y_3^6*y_4^3+11436*y_3^
       5*y_4^4+3000*y_3^4*y_4^5+2070*y_3^3*y_4^6-216*y_3^2*y_4^7-138*y_3*y_4^8-10*y_4^9-21600*y_3^8-
       11264*y_3^7*y_4+39528*y_3^6*y_4^2-38752*y_3^5*y_4^3-6144*y_3^4*y_4^4-26700*y_3^3*y_4^5-6516*y_3
       ^2*y_4^6+808*y_3*y_4^7+84*y_4^8+59184*y_3^7+147028*y_3^6*y_4+5448*y_3^5*y_4^2-126648*y_3^4*y_4^
       3+58500*y_3^3*y_4^4-30730*y_3^2*y_4^5+24336*y_3*y_4^6+4098*y_4^7+234720*y_3^6-853520*y_3^5*y_4+
       579552*y_3^4*y_4^2-146576*y_3^3*y_4^3+617328*y_3^2*y_4^4-14852*y_3*y_4^5+1812*y_4^6-825984*y_3^
       5-75220*y_3^4*y_4+452700*y_3^3*y_4^2-1441768*y_3^2*y_4^3-267768*y_3*y_4^4-249384*y_4^5+1592352*
       y_3^4+810080*y_3^3*y_4+445128*y_3^2*y_4^2-233168*y_3*y_4^3+147252*y_4^4-1449072*y_3^3+1877030*y
       _3^2*y_4+503742*y_3*y_4^2+4161430*y_4^3-4127328*y_3^2-602712*y_3*y_4-12735684*y_4^2+4200984*y_3
       +12364002*y_4-699840)







--Degenerate case 1: 

Y = matrix{{y_1},{y_2},{y_3},{y_1+y_2-y_3}}

l0 = trace(Minvdet * M0) * deter - transpose(Y) * Minvdet * M0 * Minvdet * Y
l1 = trace(Minvdet * M1) * deter - transpose(Y) * Minvdet * M1 * Minvdet * Y
l2 = trace(Minvdet * M2) * deter - transpose(Y) * Minvdet * M2 * Minvdet * Y
l3 = trace(Minvdet * M3) * deter - transpose(Y) * Minvdet * M3 * Minvdet * Y
f0 = l0_0_0
f1 = l1_0_0
f2 = l2_0_0
f3 = l3_0_0
f0R = numerator(f0)
f1R = numerator(f1)
f2R = numerator(f2)
f3R = numerator(f3)

E1 = 2*y_1^3+3*y_1^2*y_2-3*y_1^2*y_3+18*y_1*y_2*y_3-9*y_2^3+21*y_2^2*y_3-8*y_3^3

I = ideal(f0R,f1R,f2R,f3R,E1)

detR = numerator(deter)
Mn = minimalPrimes (ideal(detR))
fac1 = Mn_0_0
fac2 = Mn_1_0
J1 = saturate(I, fac1)
J2 = saturate(J1, fac2)
J2




-- s = sqrt(33), t = sqrt(9y2^2+2y2y3+9y3^2)  --> ambos adjuntados via relaciones
Rgen = QQ[y2,y3,s,t,a0,a1,a2,a3] / ideal(s^2-33, t^2-(9*y2^2+2*y2*y3+9*y3^2))

M = matrix{{a0,a1,a2,a3},
           {a1,a0,a1,a2},
           {a2,a1,a0,a1},
           {a3,a2,a1,a0}}
M0 = matrix{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}
M1 = matrix{{0,1,0,0},{1,0,1,0},{0,1,0,1},{0,0,1,0}}
M2m = matrix{{0,0,1,0},{0,0,0,1},{1,0,0,0},{0,1,0,0}}
M3 = matrix{{0,0,0,1},{0,0,0,0},{0,0,0,0},{1,0,0,0}}
deter = det M

-- adjugate manual (evita frac()/inverse() sobre anillo cociente)
adj = mutableMatrix(Rgen,4,4)
for i from 0 to 3 do (
  for j from 0 to 3 do (
    rowsToKeep := select(0..3, k -> k != j);
    colsToKeep := select(0..3, k -> k != i);
    sub4 := submatrix(M, rowsToKeep, colsToKeep);
    adj_(i,j) = (-1)^(i+j) * det(sub4);
  );
);
Minvdet = matrix adj

-- y1 vía la fórmula general del teorema para D2=0 (rama +)
y1 = (11*y3 + 3*s*t)/11
Y = matrix{{y1},{y2},{y3},{y1+y2-y3}}

l0 = trace(Minvdet*M0)*deter - transpose(Y)*Minvdet*M0*Minvdet*Y
l1 = trace(Minvdet*M1)*deter - transpose(Y)*Minvdet*M1*Minvdet*Y
l2 = trace(Minvdet*M2m)*deter - transpose(Y)*Minvdet*M2m*Minvdet*Y
l3 = trace(Minvdet*M3)*deter - transpose(Y)*Minvdet*M3*Minvdet*Y

f0 = l0_0_0
f1 = l1_0_0
f2 = l2_0_0
f3 = l3_0_0

I = ideal(f0,f1,f2,f3)
Mn = minimalPrimes ideal(deter)
fac1 = Mn_0_0
fac2 = Mn_1_0
J1 = saturate(I, fac1)
J2 = saturate(J1, fac2)


dim J2
degree J2









R = QQ[a_0, a_1, a_2, a_3,y_1, y_2, y_3, y_4]

--R = (frac QQ[y_1, y_2, y_3, y_4])[a_0, a_1, a_2, a_3]

K = frac(R)


M = matrix{{a_0, a_1, a_2, a_3},
           {a_1, a_0, a_1, a_2},
           {a_2, a_1, a_0, a_1},
           {a_3, a_2, a_1, a_0}}


M0 = matrix{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}
M1 = matrix{{0,1,0,0},{1,0,1,0},{0,1,0,1},{0,0,1,0}}
M2 = matrix{{0,0,1,0},{0,0,0,1},{1,0,0,0},{0,1,0,0}}
M3 = matrix{{0,0,0,1},{0,0,0,0},{0,0,0,0},{1,0,0,0}}

deter = det(M)


Minv = inverse(M)


Minvdet = Minv * deter


Y = matrix{{y_1},{y_2},{y_3},{y_1+y_2-y_3}}

--- 2,3,5,0
--Degenerate case 0
--Y = matrix{{-2*y_2-y_3},{y_2},{y_3},{-y_2-2*y_3}}

--Y = matrix{{-y_2},{y_2},{-y_2},{y_2}}


l0 = trace(Minvdet * M0) * deter - transpose(Y) * Minvdet * M0 * Minvdet * Y
l1 = trace(Minvdet * M1) * deter - transpose(Y) * Minvdet * M1 * Minvdet * Y
l2 = trace(Minvdet * M2) * deter - transpose(Y) * Minvdet * M2 * Minvdet * Y
l3 = trace(Minvdet * M3) * deter - transpose(Y) * Minvdet * M3 * Minvdet * Y


f0 = l0_0_0
f1 = l1_0_0
f2 = l2_0_0
f3 = l3_0_0

f0R = numerator(f0)
f1R = numerator(f1)
f2R = numerator(f2)
f3R = numerator(f3)

ring f0R
ring f1R
ring f2R
ring f3R

---E1 = 2*y_1^3+3*y_1^2*y_2-3*y_1^2*y_3+18*y_1*y_2*y_3-9*y_2^3+21*y_2^2*y_3-8*y_3^3
--FE1=numerator(E1)
--I=ideal(f0R,f1R,f2R,f3R,FE1)

I=ideal(f0R,f1R,f2R,f3R)

detR=numerator(deter)


Mn=minimalPrimes (ideal(detR))

fac1=Mn_0_0
fac2=Mn_1_0



J1 = saturate(I, fac1)


J2 = saturate(J1, fac2)


E1 = 2*y_1^3+3*y_1^2*y_2-3*y_1^2*y_3+18*y_1*y_2*y_3-9*y_2^3+21*y_2^2*y_3-8*y_3^3
FE1=numerator(E1)

newI=ideal(J2,FE1)


nJ1 = saturate(newI, fac1)


nJ2 = saturate(nJ1, fac2)

Sgb=QQ[y_1,y_2,y_3,y_4][a_3,a_2,a_1,a_0, MonomialOrder => Lex]

RJ2= substitute(nJ2, Sgb)

GRJ2 = gens gb RJ2

GRJ2_0





---Subcase extra



R = QQ[a_0, a_1, a_2, a_3,y_1, y_2, y_3, y_4]

--R = (frac QQ[y_1, y_2, y_3, y_4])[a_0, a_1, a_2, a_3]

K = frac(R)


M = matrix{{a_0, a_1, a_2, a_3},
           {a_1, a_0, a_1, a_2},
           {a_2, a_1, a_0, a_1},
           {a_3, a_2, a_1, a_0}}


M0 = matrix{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}
M1 = matrix{{0,1,0,0},{1,0,1,0},{0,1,0,1},{0,0,1,0}}
M2 = matrix{{0,0,1,0},{0,0,0,1},{1,0,0,0},{0,1,0,0}}
M3 = matrix{{0,0,0,1},{0,0,0,0},{0,0,0,0},{1,0,0,0}}

deter = det(M)


Minv = inverse(M)


Minvdet = Minv * deter


Y = matrix{{y_1},{y_2},{-y_2},{y_1+2*y_2}}

--- 2,3,5,0
--Degenerate case 0
--Y = matrix{{-2*y_2-y_3},{y_2},{y_3},{-y_2-2*y_3}}

--Y = matrix{{-y_2},{y_2},{-y_2},{y_2}}


l0 = trace(Minvdet * M0) * deter - transpose(Y) * Minvdet * M0 * Minvdet * Y
l1 = trace(Minvdet * M1) * deter - transpose(Y) * Minvdet * M1 * Minvdet * Y
l2 = trace(Minvdet * M2) * deter - transpose(Y) * Minvdet * M2 * Minvdet * Y
l3 = trace(Minvdet * M3) * deter - transpose(Y) * Minvdet * M3 * Minvdet * Y


f0 = l0_0_0
f1 = l1_0_0
f2 = l2_0_0
f3 = l3_0_0

f0R = numerator(f0)
f1R = numerator(f1)
f2R = numerator(f2)
f3R = numerator(f3)

ring f0R
ring f1R
ring f2R
ring f3R

---E1 = 2*y_1^3+3*y_1^2*y_2-3*y_1^2*y_3+18*y_1*y_2*y_3-9*y_2^3+21*y_2^2*y_3-8*y_3^3
--FE1=numerator(E1)
--I=ideal(f0R,f1R,f2R,f3R,FE1)

I=ideal(f0R,f1R,f2R,f3R)

detR=numerator(deter)


Mn=minimalPrimes (ideal(detR))

fac1=Mn_0_0
fac2=Mn_1_0



J1 = saturate(I, fac1)


J2 = saturate(J1, fac2)


E1 = 2*y_1^3+3*y_1^2*y_2-3*y_1^2*y_3+18*y_1*y_2*y_3-9*y_2^3+21*y_2^2*y_3-8*y_3^3
FE1=numerator(E1)

newI=ideal(J2,FE1)


nJ1 = saturate(newI, fac1)


nJ2 = saturate(nJ1, fac2)

Sgb=QQ[y_1,y_2,y_3,y_4][a_3,a_2,a_1,a_0, MonomialOrder => Lex]

RJ2= substitute(nJ2, Sgb)

GRJ2 = gens gb RJ2

GRJ2_0





-- case D2=11y_1^2-27y_2^2-22y_1y_3-6y_2y_3-16y_3^2=0



R = QQ[a_0, a_1, a_2, a_3,y_1, y_2, y_3, y_4]

--R = (frac QQ[y_1, y_2, y_3, y_4])[a_0, a_1, a_2, a_3]

K = frac(R)


M = matrix{{a_0, a_1, a_2, a_3},
           {a_1, a_0, a_1, a_2},
           {a_2, a_1, a_0, a_1},
           {a_3, a_2, a_1, a_0}}


M0 = matrix{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}
M1 = matrix{{0,1,0,0},{1,0,1,0},{0,1,0,1},{0,0,1,0}}
M2 = matrix{{0,0,1,0},{0,0,0,1},{1,0,0,0},{0,1,0,0}}
M3 = matrix{{0,0,0,1},{0,0,0,0},{0,0,0,0},{1,0,0,0}}

deter = det(M)


Minv = inverse(M)


Minvdet = Minv * deter


Y = matrix{{y_1},{y_2},{y_3},{y_1+y_2-y_3}}

--- 2,3,5,0
--Degenerate case 0
--Y = matrix{{-2*y_2-y_3},{y_2},{y_3},{-y_2-2*y_3}}

--Y = matrix{{-y_2},{y_2},{-y_2},{y_2}}


l0 = trace(Minvdet * M0) * deter - transpose(Y) * Minvdet * M0 * Minvdet * Y
l1 = trace(Minvdet * M1) * deter - transpose(Y) * Minvdet * M1 * Minvdet * Y
l2 = trace(Minvdet * M2) * deter - transpose(Y) * Minvdet * M2 * Minvdet * Y
l3 = trace(Minvdet * M3) * deter - transpose(Y) * Minvdet * M3 * Minvdet * Y


f0 = l0_0_0
f1 = l1_0_0
f2 = l2_0_0
f3 = l3_0_0

f0R = numerator(f0)
f1R = numerator(f1)
f2R = numerator(f2)
f3R = numerator(f3)

ring f0R
ring f1R
ring f2R
ring f3R

---E1 = 2*y_1^3+3*y_1^2*y_2-3*y_1^2*y_3+18*y_1*y_2*y_3-9*y_2^3+21*y_2^2*y_3-8*y_3^3
--FE1=numerator(E1)
--I=ideal(f0R,f1R,f2R,f3R,FE1)

I=ideal(f0R,f1R,f2R,f3R)

detR=numerator(deter)


Mn=minimalPrimes (ideal(detR))

fac1=Mn_0_0
fac2=Mn_1_0



J1 = saturate(I, fac1)


J2 = saturate(J1, fac2)


D2=11y_1^2-27y_2^2-22y_1y_3-6y_2y_3-16y_3^2=0
FD2=numerator(D2)

newI=ideal(J2,FD2)


nJ1 = saturate(newI, fac1)


nJ2 = saturate(nJ1, fac2)

Sgb=QQ[y_1,y_2,y_3,y_4][a_3,a_2,a_1,a_0, MonomialOrder => Lex]

RJ2= substitute(nJ2, Sgb)

GRJ2 = gens gb RJ2

GRJ2_0



