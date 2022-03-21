using LinearAlgebra

function Cn(v, n::Integer)
    θ = 2*π / n
    return rotation_matrix(v, θ)
end

"""
    Molecules.rotation_matrix(v, θ::AbstractFloat)

Returns a rotation matrix by θ about a rotation axis v
"""
function rotation_matrix(V, θ::AbstractFloat)
    cosθ = cos(θ)
    sinθ = sin(θ)
    v = normalize(V)
    a = [1,2,3]
    O = zeros(eltype(v), (3,3))
    O .+= 1 - cosθ
    for i = 1:3, j = 1:3
        if i == j
            O[i,i] *= v[i]^2
            O[i,i] += cosθ
        else
            O[i,j] *= v[i]*v[j]
            b = [i,j]
            C = [i for i in a if i ∉ b][1]
            if i < j
                O[i,j] += (-1)^(i+j) * v[C]*sinθ
            else
                O[i,j] += (-1)^(i+j-1) * v[C]*sinθ
            end
        end
    end
    return O
end

function rotation_matrixd(v, θ::AbstractFloat)
    r = deg2rad(θ)
    return rotation_matrix(v, r)
end

"""
    Molecules.reflection_matrix(v::Vector)

Returns a reflection matrix through the plane with normal vector v
"""
function reflection_matrix(V)
    O = zeros(eltype(V), (3,3))
    v = normalize(V)
    for i = 1:3, j = i:3
        if i == j
            O[i,i] = 1 - 2*v[i]^2
        else
            O[i,j] = -2 * v[i] * v[j]
            O[j,i] = O[i,j]
        end
    end
    return O
end

function σ(v)
    return reflection_matrix(v)
end

"""
    Molecules.Sn(v::Vector, n::int)

Returns Sn improper rotation about vector v
"""    
function Sn(v, n)
    return Cn(v, n) * σ(v)
end

"""
    Molecules.inversion_matrix()

Returns a diagonal matrix with -1 along the diagonal
"""
function inversion_matrix()
    a = zeros((3,3))
    for i = 1:3
        a[i,i] = -1
    end    
    return a
end

function i()
    return inversion_matrix()
end
    
