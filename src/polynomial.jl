#############################################################################
#############################################################################
#
# This file defines the polynomial type with several operations 
#                                                                               
#############################################################################
#############################################################################

####################################
# PolynomialDense type and construction #
####################################

"""
A PolynomialDense type - designed to be for polynomials with integer coefficients.
"""
struct PolynomialDense

    #A zero packed vector of terms
    #Terms are assumed to be in order with first term having degree 0, second degree 1, and so fourth
    #until the degree of the polynomial. The leading term (i.e. last) is assumed to be non-zero except 
    #for the zero polynomial where the vector is of length 1.
    #Note: at positions where the coefficient is 0, the power of the term is also 0 (this is how the Term type is designed)
    terms::Vector{Term}   
    
    #Inner constructor of 0 polynomial
    PolynomialDense() = new([zero(Term)])

    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialDense(vt::Vector{Term})

        #Filter the vector so that there is not more than a single zero term
        vt = filter((t)->!iszero(t), vt)
        if isempty(vt)
            vt = [zero(Term)]
        end

        max_degree = maximum((t)->t.degree, vt)
        terms = [zero(Term) for i in 0:max_degree] #First set all terms with zeros

        #now update based on the input terms
        for t in vt
            terms[t.degree + 1] = t #+1 accounts for 1-indexing
        end
        return new(terms)
    end
end

struct PolynomialSparse
    terms::Vector{Term}       
    PolynomialSparse() = new([zero(Term)])
    function PolynomialSparse(vt::Vector{Term})
        vt = filter((t)->!iszero(t), vt)
        if isempty(vt)
            vt = [zero(Term)]
        end

        terms = [] 

        for i in 1:length(vt)
            push!(terms,vt[i]) 
        end
        return new(terms)
    end
end

struct PolynomialSparseBI
    terms::Vector{TermBI}       
    PolynomialSparseBI() = new([zero(TermBI)])
    function PolynomialSparseBI(vt::Vector{TermBI})
        vt = filter((t)->!iszero(t), vt)
        if isempty(vt)
            vt = [zero(TermBI)]
        end

        terms = [] 

        for i in 1:length(vt)
            push!(terms,vt[i]) 
        end
        return new(terms)
    end
end

"""
This function maintains the invariant of the PolynomialDense type so that there are no zero terms beyond the highest
non-zero term.
"""
function trim!(p::PolynomialDense)::PolynomialDense
    i = length(p.terms)
    while i > 1
        if iszero(p.terms[i])
            pop!(p.terms)
        else
            break
        end
        i -= 1
    end
    return p
end

function trim!(p::PolynomialSparse)::PolynomialSparse
    i = length(p.terms)
    while i > 1
        if iszero(p.terms[i])
            pop!(p.terms)
        else
            break
        end
        i -= 1
    end
    return p
end

function trim!(p::PolynomialSparseBI)::PolynomialSparseBI
    i = length(p.terms)
    while i > 1
        if iszero(p.terms[i])
            pop!(p.terms)
        else
            break
        end
        i -= 1
    end
    return p
end

"""
Construct a polynomial with a single term.
"""
PolynomialDense(t::Term) = PolynomialDense([t])
PolynomialSparse(t::Term) = PolynomialSparse([t])
PolynomialSparseBI(t::TermBI) = PolynomialSparseBI([t])


"""
Construct a polynomial of the form x^p-x.
"""
cyclotonic_polynomial(p::Int) = PolynomialDense([Term(1,p), Term(-1,0)])


"""
Construct a polynomial of the form x-n.
"""
linear_monic_polynomial(n::Int) = PolynomialDense([Term(1,1), Term(-n,0)])

"""
Construct a polynomial of the form x.
"""
x_poly() = PolynomialDense(Term(1,1))

"""
Creates the zero polynomial.
"""
zero(::Type{PolynomialDense})::PolynomialDense = PolynomialDense()
zero(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse()
zero(::Type{PolynomialSparseBI})::PolynomialSparseBI = PolynomialSparseBI()


"""
Creates the unit polynomial.
"""
one(::Type{PolynomialDense})::PolynomialDense = PolynomialDense(one(Term))
one(p::PolynomialDense) = one(typeof(p))

one(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse(one(Term))
one(p::PolynomialSparse) = one(typeof(p))

one(::Type{PolynomialSparseBI})::PolynomialSparseBI = PolynomialSparseBI(one(TermBI))
one(p::PolynomialSparseBI) = one(typeof(p))


"""
Generates a random polynomial.
"""
function rand(::Type{PolynomialDense} ; 
                degree::Int = -1, 
                terms::Int = -1, 
                max_coeff::Int = 100, 
                mean_degree::Float64 = 5.0,
                prob_term::Float64  = 0.7,
                monic = false,
                condition = (p)->true)
        
    while true 
        _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
        _terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
        degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
        coeffs = rand(1:max_coeff,_terms+1)
        monic && (coeffs[end] = 1)
        p = PolynomialDense( [Term(coeffs[i],degrees[i]) for i in 1:length(degrees)] )
        condition(p) && return p
    end
end

function rand(::Type{PolynomialSparse} ; 
    degree::Int = -1, 
    terms::Int = -1, 
    max_coeff::Int = 100, 
    mean_degree::Float64 = 5.0,
    prob_term::Float64  = 0.7,
    monic = false,
    condition = (p)->true)

while true 
_degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
_terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
coeffs = rand(1:max_coeff,_terms+1)
monic && (coeffs[end] = 1)
p = PolynomialSparse( [Term(coeffs[i],degrees[i]) for i in 1:length(degrees)] )
condition(p) && return p
end
end

function rand(::Type{PolynomialSparseBI} ; 
    degree::Int = -1, 
    terms::Int = -1, 
    max_coeff::Int = 100, 
    mean_degree::Float64 = 5.0,
    prob_term::Float64  = 0.7,
    monic = false,
    condition = (p)->true)

while true 
_degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
_terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
coeffs = rand(1:max_coeff,_terms+1)
monic && (coeffs[end] = 1)
p = PolynomialSparseBI( [TermBI(coeffs[i],degrees[i]) for i in 1:length(degrees)] )
condition(p) && return p
end
end

###########
# Display #
###########

"""
Show a polynomial.
"""
function show(io::IO, p::PolynomialDense) 
    if iszero(p)
        print(io,"0")
    else
        positive_terms = []
        for (i,t) in enumerate(reverse(p.terms))
            t.coeff !=0 ? push!(positive_terms, t) : continue
        end
        for (i,t) in enumerate((positive_terms))
            if t.coeff > 0
                if t.degree == 0 
                    print(io, i != 1 ? " + " : "", t.coeff)
                end
                if t.degree == 1
                    print(io, i != 1 ? " + " : "", t.coeff, "⋅x")
                end
                if t.degree > 1
                    print(io, i != 1 ? " + " : "", t)
                end
            end
           
            if t.coeff < 0
                if t.degree ==0
                    print(io, " - ", abs(t.coeff))
                end
                if t.degree ==1
                    print(io, " - ", abs(t.coeff), "⋅x")
                end
                if t.degree > 1
                    print(io, " - ", abs(t.coeff), "⋅x^",t.degree)
                end
            end
        end
    end
end

function show(io::IO, p::PolynomialSparse) 
    if iszero(p)
        print(io,"0")
    else
        positive_terms = []
        for (i,t) in enumerate(reverse(p.terms))
            t.coeff !=0 ? push!(positive_terms, t) : continue
        end
        for (i,t) in enumerate((positive_terms))
            if t.coeff > 0
                if t.degree == 0 
                    print(io, i != 1 ? " + " : "", t.coeff)
                end
                if t.degree == 1
                    print(io, i != 1 ? " + " : "", t.coeff, "⋅x")
                end
                if t.degree > 1
                    print(io, i != 1 ? " + " : "", t)
                end
            end
           
            if t.coeff < 0
                if t.degree ==0
                    print(io, " - ", abs(t.coeff))
                end
                if t.degree ==1
                    print(io, " - ", abs(t.coeff), "⋅x")
                end
                if t.degree > 1
                    print(io, " - ", abs(t.coeff), "⋅x^",t.degree)
                end
            end
        end
    end
end

function show(io::IO, p::PolynomialSparseBI) 
    if iszero(p)
        print(io,"0")
    else
        positive_terms = []
        for (i,t) in enumerate(reverse(p.terms))
            t.coeff !=0 ? push!(positive_terms, t) : continue
        end
        for (i,t) in enumerate((positive_terms))
            if t.coeff > 0
                if t.degree == 0 
                    print(io, i != 1 ? " + " : "", t.coeff)
                end
                if t.degree == 1
                    print(io, i != 1 ? " + " : "", t.coeff, "⋅x")
                end
                if t.degree > 1
                    print(io, i != 1 ? " + " : "", t)
                end
            end
           
            if t.coeff < 0
                if t.degree ==0
                    print(io, " - ", abs(t.coeff))
                end
                if t.degree ==1
                    print(io, " - ", abs(t.coeff), "⋅x")
                end
                if t.degree > 1
                    print(io, " - ", abs(t.coeff), "⋅x^",t.degree)
                end
            end
        end
    end
end

##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the non-zero terms of the polynomial. This implements the iteration interface.
"""
iterate(p::PolynomialDense, state=1) = iterate(p.terms, state)
iterate(p::PolynomialSparse, state=1) = iterate(p.terms, state)
iterate(p::PolynomialSparseBI, state=1) = iterate(p.terms, state)


##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.
"""
length(p::Union{PolynomialDense,PolynomialSparse,PolynomialSparseBI}) = length(p.terms) 

"""
The leading term of the polynomial.
"""
leading(p::Union{PolynomialDense,PolynomialSparse})::Term = isempty(p.terms) ? zero(Term) : last(p.terms) 
leading(p::PolynomialSparseBI)::TermBI = isempty(p.terms) ? zero(TermBI) : last(p.terms) 


"""
Returns the coefficients of the polynomial.
"""
coeffs(p::Union{PolynomialDense,PolynomialSparse,PolynomialSparseBI})::Vector{Int} = [t.coeff for t in p]

"""
The degree of the polynomial.
"""
degree(p::Union{PolynomialDense,PolynomialSparse})::Int = leading(p).degree 
degree(p::PolynomialSparseBI)::Union{Int,BigInt} = leading(p).degree 


"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::Union{PolynomialDense,PolynomialSparse})::Int = euclid_alg(coeffs(p))
content(p::PolynomialSparseBI)::Union{Int,BigInt} = euclid_alg(coeffs(p))

"""
Evaluate the polynomial at a point `x`.
"""
evaluate(f::Union{PolynomialDense,PolynomialSparse,PolynomialSparseBI}, x::T) where T <: Number = sum(evaluate(t,x) for t in f)

################################
# Pushing and popping of terms #
################################

"""
Push a new term into the polynomial.
"""
#Note that ideally this would throw and error if pushing another term of degree that is already in the polynomial
function push!(p::Union{PolynomialDense,PolynomialSparse,PolynomialSparseBI}, t::Term) 
    if t.degree <= degree(p)
        p.terms[t.degree + 1] = t
    else
        append!(p.terms, zeros(Term, t.degree - degree(p)-1))
        push!(p.terms, t)
    end
    return p        
end

"""
Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.
"""
function pop!(p::Union{PolynomialDense,PolynomialSparse,PolynomialSparseBI})::Term 
    popped_term = pop!(p.terms) #last element popped is leading coefficient

    while !isempty(p.terms) && iszero(last(p.terms))
        pop!(p.terms)
    end

    if isempty(p.terms)
        push!(p.terms, zero(Term))
    end

    return popped_term
end

"""
Check if the polynomial is zero.
"""
iszero(p::Union{PolynomialDense,PolynomialSparse,PolynomialSparseBI})::Bool = p.terms == [Term(0,0)]

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

"""
The negative of a polynomial.
"""
-(p::PolynomialDense) = PolynomialDense(map((pt)->-pt, p.terms))
-(p::PolynomialSparse) = PolynomialSparse(map((pt)->-pt, p.terms))
-(p::PolynomialSparseBI) = PolynomialSparseBI(map((pt)->-pt, p.terms))


"""
Create a new polynomial which is the derivative of the polynomial.
"""
function derivative(p::PolynomialDense)::PolynomialDense 
    der_p = PolynomialDense()
    for term in p
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p,der_term)
    end
    return trim!(der_p)
end

function derivative(p::PolynomialSparse)::PolynomialSparse 
    der_p = PolynomialSparse()
    for term in p
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p,der_term)
    end
    return trim!(der_p)
end

function derivative(p::PolynomialSparseBI)::PolynomialSparseBI 
    der_p = PolynomialSparseBI()
    for term in p
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p,der_term)
    end
    return trim!(der_p)
end

"""
The prim part (multiply a polynomial by the inverse of its content).
"""
prim_part(p::Union{PolynomialDense,PolynomialSparse,PolynomialSparseBI}) = p ÷ content(p)


"""
A square free polynomial.
"""
square_free(p::PolynomialDense, prime::Int)::PolynomialDense = (p ÷ gcd(p,derivative(p),prime))(prime)
square_free(p::PolynomialSparse, prime::Int)::PolynomialSparse = (p ÷ gcd(p,derivative(p),prime))(prime)
square_free(p::PolynomialSparseBI, prime::BigInt)::PolynomialSparseBI = (p ÷ gcd(p,derivative(p),prime))(prime)


#################################
# Queries about two polynomials #
#################################

"""
Check if two polynomials are the same
"""
==(p1::PolynomialDense, p2::PolynomialDense)::Bool = p1.terms == p2.terms
==(p1::PolynomialSparse, p2::PolynomialSparse)::Bool = p1.terms == p2.terms
==(p1::PolynomialSparseBI, p2::PolynomialSparseBI)::Bool = p1.terms == p2.terms

"""
Check if a polynomial is equal to 0.
"""
#Note that in principle there is a problem here. E.g The polynomial 3 will return true to equalling the integer 2.
==(p::PolynomialDense, n::T) where T <: Real = iszero(p) == iszero(n)
==(p::PolynomialSparse, n::T) where T <: Real = iszero(p) == iszero(n)
==(p::PolynomialSparseBI, n::T) where T <: Real = iszero(p) == iszero(n)


##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

"""
Subtraction of two polynomials.
"""
-(p1::PolynomialDense, p2::PolynomialDense)::PolynomialDense = p1 + (-p2)
-(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse = p1 + (-p2)
-(p1::PolynomialSparseBI, p2::PolynomialSparseBI)::PolynomialSparse = p1 + (-p2)


"""
Multiplication of polynomial and term.
"""
*(t::Term, p1::PolynomialDense)::PolynomialDense = iszero(t) ? PolynomialDense() : PolynomialDense(map((pt)->t*pt, p1.terms))
*(p1::PolynomialDense, t::Term)::PolynomialDense = t*p1

*(t::Term, p1::PolynomialSparse)::PolynomialSparse = iszero(t) ? PolynomialSparse() : PolynomialSparse(map((pt)->t*pt, p1.terms))
*(p1::PolynomialSparse, t::Term)::PolynomialSparse = t*p1

*(t::TermBI, p1::PolynomialSparseBI)::PolynomialSparseBI = iszero(t) ? PolynomialSparseBI() : PolynomialSparseBI(map((pt)->t*pt, p1.terms))
*(p1::PolynomialSparseBI, t::TermBI)::PolynomialSparseBI = t*p1


"""
Multiplication of polynomial and an integer.
"""
*(n::Int, p::PolynomialDense)::PolynomialDense = p*Term(n,0)
*(p::PolynomialDense, n::Int)::PolynomialDense = n*p

*(n::Int, p::PolynomialSparse)::PolynomialSparse = p*Term(n,0)
*(p::PolynomialSparse, n::Int)::PolynomialSparse = n*p

*(n::BigInt, p::PolynomialSparseBI)::PolynomialSparseBI = p*TermBI(n,0)
*(p::PolynomialSparseBI, n::BigInt)::PolynomialSparseBI = n*p

"""
Integer division of a polynomial by an integer.

Warning this may not make sense if n does not divide all the coefficients of p.
"""
÷(p::PolynomialDense, n::Int) = (prime)->PolynomialDense(map((pt)->((pt ÷ n)(prime)), p.terms))
÷(p::PolynomialSparse, n::Int) = (prime)->PolynomialSparse(map((pt)->((pt ÷ n)(prime)), p.terms))
÷(p::PolynomialSparseBI, n::BigInt) = (prime)->PolynomialSparseBI(map((pt)->((pt ÷ n)(prime)), p.terms))


"""
Take the mod of a polynomial with an integer.
"""
function mod(f::Union{PolynomialDense,PolynomialSparse}, p::Int)::Union{PolynomialDense,PolynomialSparse}
    f_out = deepcopy(f)
    for i in 1:length(f_out.terms)
        f_out.terms[i] = mod(f_out.terms[i], p)
    end
    return trim!(f_out)
end

function mod(f::PolynomialSparseBI, p::BigInt)::PolynomialSparseBI
    f_out = deepcopy(f)
    for i in 1:length(f_out.terms)
        f_out.terms[i] = mod(f_out.terms[i], p)
    end
    return trim!(f_out)
end

"""
Power of a polynomial mod prime.
"""
function pow_mod(p::Union{PolynomialDense,PolynomialSparse}, n::Int, prime::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
        out = mod(out, prime)
    end
    return out
end

function pow_mod(p::PolynomialSparseBI, n::BigInt, prime::BigInt)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
        out = mod(out, prime)
    end
    return out
end