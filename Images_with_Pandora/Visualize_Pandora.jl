using Pkg
using Oscar

#Pkg.add(url="https://github.com/tbrysiewicz/Pandora.git")

using Pandora


#Put variables using in the others code

E= EnumerativeProblem([fa1_expr,fa2_expr,fa3_expr,fa4_expr],variables=[a1,a2,a3,a4],parameters=vcat(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10),torus_only=true)



degree(E)

bkk_bound(E)

bezout_bound(E)

knowledge(E)

last(knowledge(E))

#G = monodromy_group(E)


#is_primitive(G)

#order(G)

#gens(G)

#is_decomposable(E)

is_lacunary(E)

ER = tally.(explore(E, [n_real_solutions, n_positive_solutions]), n_samples = 100)


O = maximize_n_real_solutions(E)

C = certify(record_fibre(O),E)


(V,P) = visualize(E; near = record_parameters(O), strategy = :quadtree)

using Pandora
n_total_solutions(solutions) = length(solutions)

my_n_real_solutions(solutions) = count(x -> all(isreal, x), solutions)

my_n_positive_solutions(solutions) = count(x -> all(isreal, x) && all(y -> y > 0, x), solutions)

props = [n_total_solutions, my_n_real_solutions, my_n_positive_solutions]

results = explore(E, props; n_samples=1000)

# Conteo de resultados
ER = tally.(results)

# Mostrar resultados
println(ER)

ER = tally.(explore(E, props, n_samples = 1000))

O = maximize_n_real_solutions(E)
C = certify(record_fibre(O),E)


(V,P) = Pandora.visualize(E; near = record_parameters(O), strategy = :quadtree)

refine!(V);
refine!(V);
P = Pandora.visualize(V)



#mkpath("/home/gabriel/OutputFiles")


#Pandora.save(P, "Ar4refine.png")
