struct Symel
    symbol
    rrep
end

struct Op
    optype
    order
    nprimes
end

struct Chartable
    name
    symels
    irreps
    characters
end

struct PG
    str
    family
    n
    subfamily
end

struct Rep
    irrep::Bool
    d::Int
    p::Int
    v::Int
    h::Int
    i::Int
    l::Int
    name
end

function Base.:(==)(A::Symel, B::Symel)
    if sum(A.rrep .- B.rrep) < tol
        return true
    else
        return false
    end
end