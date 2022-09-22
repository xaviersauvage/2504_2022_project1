#############################################################################
#############################################################################
#
# This file defines the Term type with several operations 
#                                                                               
#############################################################################
#############################################################################

##############################
# Term type and construction #
##############################

"""
A term.
"""
struct Term  #structs are immutable by default
    coeff::Int
    degree::Int
    function Term(coeff::Int, degree::Int)
        degree < 0 && error("Degree must be non-negative")
        coeff != 0 ? new(coeff,degree) : new(coeff,0)
    end
end

struct TermBI  #structs are immutable by default
    coeff::BigInt
    degree::BigInt
    function TermBI(coeff::BigInt, degree::BigInt)
        degree < 0 && error("Degree must be non-negative")
        coeff != 0 ? new(coeff,degree) : new(coeff,0)
    end
end

"""
Creates the zero term.
"""
zero(::Type{Term})::Term = Term(0,0)
zero(::Type{TermBI})::TermBI = TermBI(big(0),big(0))

"""
Creates the unit term.
"""
one(::Type{Term})::Term = Term(1,0)
one(::Type{TermBI})::TermBI = TermBI(1,0)

###########
# Display #
###########

"""
Show a term.
"""
show(io::IO, t::Term) = print(io, "$(t.coeff)⋅x^$(t.degree)") #\cdot + [TAB]
show(io::IO, t::TermBI) = print(io, "$(t.coeff)⋅x^$(t.degree)") #\cdot + [TAB]

########################
# Queries about a term #
########################

"""
Check if a term is 0.
"""
iszero(t::Term)::Bool = iszero(t.coeff)
iszero(t::TermBI)::Bool = iszero(t.coeff)


"""
Compare two terms.
"""
isless(t1::Term,t2::Term)::Bool =  t1.degree == t2.degree ? (t1.coeff < t2.coeff) : (t1.degree < t2.degree) 
isless(t1::TermBI,t2::TermBI)::Bool =  t1.degree == t2.degree ? (t1.coeff < t2.coeff) : (t1.degree < t2.degree)  


"""
Evaluate a term at a point x.
"""
evaluate(t::Term, x::T) where T <: Number = t.coeff * x^t.degree
evaluate(t::TermBI, x::T) where T <: Number = t.coeff * x^t.degree


##########################
# Operations with a term #
##########################

"""
Add two terms of the same degree.
"""
function +(t1::Term,t2::Term)::Term 
    if t1 == zero(Term) || t2 == zero(Term)
        Term(t1.coeff + t2.coeff, maximum([t1.degree, t2.degree]))
    else
    @assert t1.degree == t2.degree
    Term(t1.coeff + t2.coeff, t1.degree)
    end
end

function +(t1::TermBI,t2::TermBI)::TermBI
    if t1 == zero(TermBI) || t2 == zero(TermBI)
        TermBI(t1.coeff + t2.coeff, maximum([t1.degree, t2.degree]))
    else
    @assert t1.degree == t2.degree
    TermBI(t1.coeff + t2.coeff, t1.degree)
    end
end

"""
Negate a term.
"""
-(t::Term,) = Term(-t.coeff,t.degree)
-(t::TermBI,) = TermBI(-t.coeff,t.degree)  


"""
Subtract two terms with the same degree.
"""
-(t1::Term, t2::Term)::Term = t1 + (-t2)
-(t1::TermBI, t2::TermBI)::TermBI = t1 + (-t2) 


"""
Multiply two terms.
"""
*(t1::Term, t2::Term)::Term = Term(t1.coeff * t2.coeff, t1.degree + t2.degree)
*(t1::TermBI, t2::TermBI)::TermBI = TermBI(t1.coeff * t2.coeff, t1.degree + t2.degree)



"""
Compute the symmetric mod of a term with an integer.
"""
mod(t::Term, p::Int) = Term(mod(t.coeff,p), t.degree)
mod(t::TermBI, p::BigInt) = TermBI(mod(t.coeff,p), t.degree)


"""
Compute the derivative of a term.
"""
derivative(t::Term) = Term(t.coeff*t.degree,max(t.degree-1,0))
derivative(t::TermBI) = TermBI(t.coeff*t.degree,max(t.degree-1,0))

"""
Divide two terms. Returns a function of an integer.
"""
function ÷(t1::Term,t2::Term) #\div + [TAB]
    @assert t1.degree ≥ t2.degree
    f(p::Int)::Term = Term(mod((t1.coeff * int_inverse_mod(t2.coeff, p)), p), t1.degree - t2.degree)
end

function ÷(t1::TermBI,t2::TermBI) #\div + [TAB]
    @assert t1.degree ≥ t2.degree
    f(p::BigInt)::TermBI = TermBI(mod((t1.coeff * int_inverse_mod(t2.coeff, p)), p), t1.degree - t2.degree)
end

"""
Integer divide a term by an integer.
"""
÷(t::Term, n::Int) = t ÷ Term(n,0)
÷(t::TermBI, n::BigInt) = t ÷ TermBI(n,0)

"""
Check if two terms are the same
"""
==(t1::Term, t2::Term)::Bool = t1.coeff == t2.coeff && t1.degree == t2.degree
==(t1::TermBI, t2::TermBI)::Bool = t1.coeff == t2.coeff && t1.degree == t2.degree
