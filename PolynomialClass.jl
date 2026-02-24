#=
    The main file for computing the AGP as objects. Note that this is set up for 2 Degrees of Freedom (x and y) and one perturbation (β).
    It shouldn't be too difficult be expand this to higher dimensions or different perturbations. 
    Note that the inputs are formatted (x, px, y, py, β)
=#
import Base: *, /, +, -

struct Term
    coeff::Float64
    powers::NTuple{5, Int}
end

struct Polynomial
    terms::Vector{Term}
end

function polynomial_to_function(poly::Polynomial)
    if isempty(poly.terms)
        return (x, px, y, py, β) -> 0.0
    end
    return function (x, px, y, py, β)
        vars = @SVector [x, px, y, py, β]
        s = 0.0
        @inbounds for term in poly.terms
            prod = term.coeff
            @inbounds for j in 1:5
                p = term.powers[j]
                if p == 0
                    continue
                elseif p == 1
                    prod *= vars[j]
                else
                    prod *= vars[j]^p
                end
            end
            s += prod
        end
        return s
    end
end

function Derivative(term::Term, var::Symbol)
    index = Dict(:x=> 1, :px=>2, :y=>3, :py=>4, :β=>5)[var]
    if term.powers[index] == 0
        return nothing
    end
    new_coeff = term.coeff * term.powers[index]

    new_powers = Base.setindex(term.powers, term.powers[index] - 1, index)
    return Term(new_coeff, new_powers)
end

function Polynomial_Derivative(poly::Polynomial, var::Symbol)
    new_terms = Term[]
    for term in poly.terms
        der_term = Derivative(term, var)
        if !isnothing(der_term)
            push!(new_terms, der_term)
        end
    end
    return simplify(Polynomial(new_terms))
end

function multiply(a::Term, b::Term)
    new_coeff = a.coeff * b.coeff
    new_powers = ntuple(i -> a.powers[i] + b.powers[i], 5)
    return Term(new_coeff, new_powers)
end

function simplify(poly::Polynomial)
    if isempty(poly.terms)
        return Polynomial([])
    end
    new_terms = Dict{NTuple{5, Int}, Float64}()
    for term in poly.terms
        key = term.powers  # already a tuple
        if haskey(new_terms, key)
            new_terms[key] += term.coeff
        else
            new_terms[key] = term.coeff
        end
    end
    simplified_terms = []
    for (powers, coeff) in new_terms
        if coeff != 0.0
            push!(simplified_terms, Term(coeff, powers))  # powers is already a tuple
        end
    end
    return Polynomial(simplified_terms)
end

function AGP_PB(poly::Polynomial)
    Ham = [
        Term(1.0, (1, 0, 0, 0, 0)),
        Term(1.0, (1, 0, 2, 0, 1)),
        Term(-1.0,(0, 1, 0, 0, 0)),
        Term(1.0, (0, 0, 1, 0, 0)),
        Term(1.0, (2, 0, 1, 0, 1)),
        Term(-1.0,(0, 0, 0, 1, 0))
    ]
    DeriVar = [:px,:px,:x,:py,:py,:y]
    new_terms = Term[]
    for term in poly.terms
        for (ham_term,var) in zip(Ham, DeriVar)
            derivative_term = Derivative(term, var)
            if !isnothing(derivative_term)
                push!(new_terms, multiply(derivative_term,ham_term))
            end
        end
    end
    return simplify(Polynomial(new_terms))
end

function nestedPB(poly::Polynomial, k::Int)
    result = poly
    for i in 1:k
        result = AGP_PB(result)
    end
    return result
end

function PB(A::Polynomial, B::Polynomial)
    vars = [:x, :px, :y, :py]
    idx = Dict(:x=>1, :px=>2, :y=>3, :py=>4)
    result_terms = Term[]
    for (v1, v2) in [(:x, :px), (:y, :py)]
        for a in A.terms
            da = Derivative(a, v1)
            if isnothing(da) continue end
            for b in B.terms
                db = Derivative(b, v2)
                if isnothing(db) continue end
                push!(result_terms, multiply(da, db))
            end
        end
        for a in A.terms
            da = Derivative(a, v2)
            if isnothing(da) continue end
            for b in B.terms
                db = Derivative(b, v1)
                if isnothing(db) continue end
                push!(result_terms, multiply(Term(-1.0, (0,0,0,0,0)), multiply(da, db)))
            end
        end
    end
    return simplify(Polynomial(result_terms))
end

function +(p1::Polynomial, p2::Polynomial)
    return simplify(Polynomial(vcat(p1.terms, p2.terms)))
end


function -(p1::Polynomial, p2::Polynomial)
    neg_terms = [Term(-term.coeff, term.powers) for term in p2.terms]
    return simplify(Polynomial(vcat(p1.terms, neg_terms)))
end

function *(a::Number, poly::Polynomial)
    new_terms = [Term(term.coeff * a, term.powers) for term in poly.terms]
    return simplify(Polynomial(new_terms))
end

function *(poly::Polynomial, a::Number)
    return a * poly
end

function /(poly::Polynomial, a::Number)
    new_terms = [Term(term.coeff / a, term.powers) for term in poly.terms]
    return simplify(Polynomial(new_terms))
end
function Print_Polynomial(poly::Polynomial)
    terms_str = String[]
    for term in poly.terms
        if term.coeff == 0.0
            continue
        end
        str = "$(term.coeff)"
        vars = ["x", "px", "y", "py", "β"]
        for (v, p) in zip(vars, term.powers)
            if p == 0
                continue
            elseif p == 1
                str *= " * $v"
            else
                str *= " * $v^$p"
            end
        end
        push!(terms_str, str)
    end
    println(join(terms_str, " + "))
end
