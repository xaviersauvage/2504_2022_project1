#############################################################################
#############################################################################
#
# This file implements polynomial addition 
#                                                                               
#############################################################################
#############################################################################

"""
Add a polynomial and a term.
"""
function +(p::PolynomialDense, t::Term)
    p = deepcopy(p)
    if t.degree > degree(p)
        push!(p, t)
    else
        if !iszero(p.terms[t.degree + 1]) #+1 is due to indexing
            p.terms[t.degree + 1] += t
        else
            p.terms[t.degree + 1] = t
        end
    end

    return trim!(p)
end

function +(p::PolynomialSparse, t::Term)
    p1 = deepcopy(p)
    t1 = deepcopy(t)
    max_degreep1 = maximum((t)->t.degree, p1.terms)
    degreet1 = t1.degree
    max_both = maximum([max_degreep1,degreet1])

    termsp1 = [zero(Term) for i  in 0:max_both]
    termst1 = [zero(Term) for i  in 0:max_both]
    termst1[t1.degree + 1] = t
    
    for t in p1.terms
        termsp1[t.degree + 1] = t
    end

    vt = termsp1 + termst1
    return PolynomialSparse(vt)

    return trim!(p)
end

function +(p::PolynomialSparseBI, t::TermBI)
    p1 = deepcopy(p)
    t1 = deepcopy(t)
    max_degreep1 = maximum((t)->t.degree, p1.terms)
    degreet1 = t1.degree
    max_both = maximum([max_degreep1,degreet1])

    termsp1 = [zero(TermBI) for i  in 0:max_both]
    termst1 = [zero(TermBI) for i  in 0:max_both]
    termst1[t1.degree + 1] = t
    
    for t in p1.terms
        termsp1[t.degree + 1] = t
    end
    vt = termsp1 + termst1
    return PolynomialSparseBI(vt)

    return trim!(p)
end

+(t::Union{Term,TermBI}, p::Union{PolynomialDense,PolynomialSparse,PolynomialSparseBI}) = p + t

"""
Add two polynomials.
"""
function +(p1::PolynomialDense, p2::PolynomialDense)::PolynomialDense
    p = deepcopy(p1)
    for t in p2
        p += t
    end
    return p
end

function +(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    p = deepcopy(p1)
    
    for t in p2
        if t.degree ∉ [t.degree for t in p]
            push!(p.terms,t)
        else
            p += t
        end
    end
    return p
end


function +(p1::PolynomialSparseBI, p2::PolynomialSparseBI)::PolynomialSparseBI
    pa = deepcopy(p1)
    pb = deepcopy(p2)
    max_degreep1 = maximum((t)->t.degree, pa.terms)
    max_degreep2 = maximum((t)->t.degree, pb.terms)
    max_both = maximum([max_degreep1,max_degreep2])

    termsp1 = [zero(TermBI) for i  in 0:max_both]
    termsp2 = [zero(TermBI) for i  in 0:max_both]

    
    for t in pa.terms
        termsp1[t.degree + 1] = t
    end

    for t in pb.terms
        termsp2[t.degree + 1] = t
    end
    vt = termsp1 + termsp2
    return PolynomialSparseBI(vt)
end

"""
Add a polynomial and an integer.
"""
+(p::PolynomialDense, n::Int) = p + Term(n,0)
+(p::PolynomialSparse, n::Int) = p + Term(n,0)
+(p::PolynomialSparseBI, n::BigInt) = p + TermBI(big(n),big(0))


+(n::Int, p::PolynomialDense) = p + Term(n,0)
+(n::Int, p::PolynomialSparse) = p + Term(n,0)
+(n::BigInt, p::PolynomialSparseBI) = p + TermBI(big(n),big(0))