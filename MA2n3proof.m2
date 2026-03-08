
R = (frac QQ[y_1, y_2, y_3])[a_0, a_1, a_2]


K = frac(R)


M = matrix{{a_0, a_1, a_2},
           {a_1, a_0, a_1},
           {a_2, a_1, a_0}}


M0 = matrix{{1,0,0},{0,1,0},{0,0,1}}
M1 = matrix{{0,1,0},{1,0,1},{0,1,0}}
M2 = matrix{{0,0,1},{0,0,0},{1,0,0}}


deter = det(M)


Minv = inverse(M)


Minvdet = Minv * deter


Y = matrix{{y_1},{y_2},{y_3}}


l0 = trace(Minvdet * M0) * deter - transpose(Y) * Minvdet * M0 * Minvdet * Y
l1 = trace(Minvdet * M1) * deter - transpose(Y) * Minvdet * M1 * Minvdet * Y
l2 = trace(Minvdet * M2) * deter - transpose(Y) * Minvdet * M2 * Minvdet * Y


f0 = l0_0_0
f1 = l1_0_0
f2 = l2_0_0


f0R = numerator(f0)
f1R = numerator(f1)
f2R = numerator(f2)


ring f0R
ring f1R
ring f2R

I=ideal(f0R,f1R,f2R)

detR=numerator(deter)

Mn=minimalPrimes (ideal(detR))

fac1=Mn_0_0
fac2=Mn_1_0



J1 = saturate(I, fac1)


J2 = saturate(J1, fac2)


RJ=radical(J2)

RJ2=RJ

MPRJ= minimalPrimes RJ


S= (frac QQ[y_1, y_2, y_3])[a_2, a_1, a_0, MonomialOrder => Lex]

RJ2= substitute(RJ, S)
GRJ2 = gens gb RJ2

GRJ2_0
GRJ2_1
GRJ2_4


S= (frac QQ[y_1, y_2, y_3])[a_0, a_1, a_2, MonomialOrder => Lex]

RJ2= substitute(RJ, S)
GRJ2 = gens gb RJ2


S= (frac QQ[y_1, y_2, y_3])[a_1, a_2, a_0, MonomialOrder => Lex]

RJ2= substitute(RJ, S)
GRJ2 = gens gb RJ2

--(4y_1^2-16y_2^2+8y_1y_3+4y_3^2)a_2^2+(-4y_1^2y_2^2+8y_2^4-4y_1^3y_3+16y_1y_2^2y_3-8y_1^2y_3^2-4y_2^2y_3^2-4y_1y_3^3)a_2-3y_1^2y_2^4+2y_1^3y_2^2y_3+2y_1y_2^4y_3+y_1^4y_3^2-4y_1^2y_2^2y_3^2-3y_2^4y_3^2+2y_1^3y_3^3+2y_1y_2^2y_3^3+y_1^2y_3^4 |



--| 4y_2a_1+(-2y_1-2y_3)a_2-y_1y_2^2+y_1^2y_3-y_2^2y_3+y_1y_3^2 |

--| 3y_2a_0+(-y_1-y_3)a_1-y_2a_2-y_2^3+y_1y_2y_3 |

--sub(GRJ2, {y_1=>3, y_2=>5, y_3=>7})

--GRJ2_0
--sub(GRJ2_0, {y_1=>2, y_2=>3, y_3=>7})
--GRJ2_1
--sub(GRJ2_1, {y_1=>2, y_2=>3, y_3=>7})
--GRJ2_4
--sub(GRJ2_4, {y_1=>2, y_2=>3, y_3=>7})



--Result A2 

| (4y_1^2-16y_2^2+8y_1y_3+4y_3^2)a_2^2+(-4y_1^2y_2^2+8y_2^4-4y_1^3y_3+16y_1y_2^2y_3-8y_1^2y_3^2-4y_2^2y_3^2-4y_1y_3^3)a_2-3y_1^2y_2^4+2y_1^3y_2^2y_3+2y_1y_2^4y_3+y_1^4y_3^2-4y_1^2y_2^2y_3^2-3y_2^4y_3^2+2y_1^3y_3^3+2y_1y_2^2y_3^3+y_1^2y_3^4 |


--Result A1
| (2y_1^2-8y_2^2+4y_1y_3+2y_3^2)a_1^2+(-2y_1^3y_2+6y_1y_2^3-4y_1^2y_2y_3+6y_2^3y_3-4y_1y_2y_3^2-2y_2y_3^3)a_1-y_1^2y_2^4+y_1^3y_2^2y_3-2y_1y_2^4y_3+2y_1^2y_2^2y_3^2-y_2^4y_3^2+y_1y_2^2y_3^3 |


--Result A0
| (12y_1^2-48y_2^2+24y_1y_3+12y_3^2)a_0^2+(-4y_1^4+40y_2^4-8y_1^3y_3-8y_1y_2^2y_3-8y_1^2y_3^2-8y_1y_3^3-4y_3^4)a_0-3y_1^2y_2^4-8y_2^6+2y_1^3y_2^2y_3+2y_1y_2^4y_3+y_1^4y_3^2+4y_1^2y_2^2y_3^2-3y_2^4y_3^2+2y_1^3y_3^3+2y_1y_2^2y_3^3+y_1^2y_3^4 |




--Result A0
(12y_1^2-48y_2^2+24y_1y_3+12y_3^2)a_0^2+(-4y_1^4+40y_2^4-8y_1^3y_3-8y_1y_2^2y_3-8y_1^2y_3^2-8y_1y_3^3-4y_3^4)a_0-3y_1^2y_2^4-8y_2^6+2y_1^3y_2^2y_3+2y_1y_2^4y_3+y_1^4y_3^2+4y_1^2y_2^2y_3^2-3y_2^4y_3^2+2y_1^3y_3^3+2y_1y_2^2y_3^3+y_1^2y_3^4 

A0=(12*y_1^2-48*y_2^2+24*y_1*y_3+12*y_3^2)
B0=(-4*y_1^4+40*y_2^4-8*y_1^3*y_3-8*y_1*y_2^2*y_3-8*y_1^2*y_3^2-8*y_1*y_3^3-4*y_3^4)
C0=-3*y_1^2*y_2^4-8*y_2^6+2*y_1^3*y_2^2*y_3+2*y_1*y_2^4*y_3+y_1^4*y_3^2+4*y_1^2*y_2^2*y_3^2-3*y_2^4*y_3^2+2*y_1^3*y_3^3+2*y_1*y_2^2*y_3^3+y_1^2*y_3^4 

discriminantA0=B0^2 - 4*A0*C0

(-B0-sqrt(discriminantA0))/(2*A0)

--Result A1
(2y_1^2-8y_2^2+4y_1y_3+2y_3^2)a_1^2+(-2y_1^3y_2+6y_1y_2^3-4y_1^2y_2y_3+6y_2^3y_3-4y_1y_2y_3^2-2y_2y_3^3)a_1-y_1^2y_2^4+y_1^3y_2^2y_3-2y_1y_2^4y_3+2y_1^2y_2^2y_3^2-y_2^4y_3^2+y_1y_2^2y_3^3 |

A1=(2*y_1^2-8*y_2^2+4*y_1*y_3+2*y_3^2)
B1=-2*y_1^3*y_2+6*y_1*y_2^3-4*y_1^2*y_2*y_3+6*y_2^3*y_3-4*y_1*y_2*y_3^2-2*y_2*y_3^3
C1=-y_1^2*y_2^4+y_1^3*y_2^2*y_3-2*y_1*y_2^4*y_3+2*y_1^2*y_2^2*y_3^2-y_2^4*y_3^2+y_1*y_2^2*y_3^3 

discriminantA1=B1^2 - 4*A1*C1
---Result A2

(4y_1^2-16y_2^2+8y_1y_3+4y_3^2)a_2^2+(-4y_1^2y_2^2+8y_2^4-4y_1^3y_3+16y_1y_2^2y_3-8y_1^2y_3^2-4y_2^2y_3^2-4y_1y_3^3)a_2-3y_1^2y_2^4+2y_1^3y_2^2y_3+2y_1y_2^4y_3+y_1^4y_3^2-4y_1^2y_2^2y_3^2-3y_2^4y_3^2+2y_1^3y_3^3+2y_1y_2^2y_3^3+y_1^2y_3^4 

A2=4*y_1^2-16*y_2^2+8*y_1*y_3+4*y_3^2
B2=-4*y_1^2*y_2^2+8*y_2^4-4*y_1^3*y_3+16*y_1*y_2^2*y_3-8*y_1^2*y_3^2-4*y_2^2*y_3^2-4*y_1*y_3^3
C2=-3*y_1^2*y_2^4+2*y_1^3*y_2^2*y_3+2*y_1*y_2^4*y_3+y_1^4*y_3^2-4*y_1^2*y_2^2*y_3^2-3*y_2^4*y_3^2+2*y_1^3*y_3^3+2*y_1*y_2^2*y_3^3+y_1^2*y_3^4 

discriminantA2=B2^2 - 4*A2*C2




S= (frac QQ[y_1, y_2, y_3])[a_2, a_1, a_0, MonomialOrder => Lex]

----
R0 = QQ[y1, y2, y3, a2, a1, a0, MonomialOrder => Lex]

ra0R = substitute(ra0, R)

I = ideal(ra0R)

Ja0 = eliminate({a0}, I)  -- eliminate a0

gens Ja0  -- polynomial a1, a2, y1, y2, y3

A2=(4*y_1^2-16*y_2^2+8*y_1y_3+4*y_3^2)
B2=(-4*y_1^2y_2^2+8*y_2^4-4*y_1^3*y_3+16*y_1y_2^2*y_3-8*y_1^2*y_3^2-4*y_2^2*y_3^2-4*y_1y_3^3)
C2=-3*y_1^2y_2^4+2*y_1^3y_2^2y_3+2*y_1y_2^4y_3+y_1^4y_3^2-4y_1^2y_2^2y_3^2-3*y_2^4y_3^2+2*y_1^3y_3^3+2*y_1y_2^2)(y_3^3)+(y_1^2)(y_3^4) 

A = coefficient(ra0, a0, 2)
B = coefficient(ra0, a0, 1)
C = coefficient(ra0, a0, 0)

