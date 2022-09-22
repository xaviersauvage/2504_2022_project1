using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")

p1 = rand(PolynomialSparse)
p2 = rand(PolynomialSparse)
p3 = rand(PolynomialSparse)
p4 = rand(PolynomialSparse)

println("p1 = ", p1)
println("p2 = ", p2)
println("p3 = ", p3)
println("p4 = ", p4)

println("\nAddition")

println("p1 + p2 = ", p1 + p2)
println("p1 + p2 + p3 + p4 = ", p1 + p2 + p3 + p4)
println("p1 + x^2 = ",p1 + Term(1,2))
println("p1 + 3 = ", p1 + 3)

println("\nMultiplication")

println("p1*p2 = ", p1 + p2)
println("p1*p2*p3*p4 = ", p1 + p2 + p3 + p4)
println("p1*x^2 = ",p1 + Term(1,2))
println("p1*3 = ", p1 + 3)
println("p1^3 = ", p1^3)

println("\nThe extended euclidean algorithm for polynomials modulo prime.")
println("Extended euclidean algorthim for p4 and p4 mod 19 = ", extended_euclid_alg(p4, p3, 19))


println("\nMod")
println("p3 mod 3 = ", mod(p3,3))
println("p4 mod 11 = ", mod(p4,11))


println("\nThe GCD of two polynomials modulo prime")
println("GCD(p1,p2) modulo 101 = ", gcd(p1,p2,101) )

println("\nDivision of polynomials modulo prime")
println("divide(p1,p2)(23) = ", divide(p1,p2)(23) )

println("\nFactorisation")
println("factor(p4,5) = ", factor(p4,5))