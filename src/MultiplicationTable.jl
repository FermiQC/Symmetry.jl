
struct MTable
    symels
    table
end

struct Group
    idxs
    elements
    mtable
    classes
    order
    ctable
end

function multifly(symels, A::Symel, B::Symel)
    Crrep = A.rrep * B.rrep
    for (i,g) in enumerate(symels)
        s = sum(abs.(g.rrep - Crrep))
        if abs(s) < tol
            return i,g
        end
    end
    println(i,g)
    throw(ArgumentError("No match found!"))
end

function build_mult_table(symels)
    h = length(symels)
    t = zeros(Int64,h,h)
    for (i,a) in enumerate(symels)
        for (j,b) in enumerate(symels)
            t[i,j] = multifly(symels, a, b)[1]
        end
    end
    return MTable(symels, t)
end

function find_inv(mtable, g)
    for i = 1:length(mtable.symels)
        if mtable.table[i,g] == 1
            return i
        end
    end
end

function find_classes(mtable)
    h = length(mtable.symels)
    assigned_elements = []
    classes = []
    for a = 1:h
        if a ∉ assigned_elements
            results = []
            for b = 1:h
                if a != b
                    binv = find_inv(mtable, b)
                    r1 = mtable.table[binv,a]
                    r = mtable.table[r1,b]
                    push!(results, r)
                end
            end
            push!(classes,unique(results))
            for i in unique(results)
                push!(assigned_elements, i)
            end
        end
    end
    return classes
end

function crst(mtable, R, S, T)
    t = T[1] # Choice of t doesn't matter
    counter = 0
    for r in R
        for s in S
            if mtable.table[r,s] == t
                counter += 1
            end
        end
    end
    return counter
end

function build_Mk(mtable, classes)
    N = length(classes)
    M = []
    for k = 1:N
        Mk = zeros(Float64,N,N)
        for i = 1:N
            for j = 1:N
                Mk[i,j] = crst(mtable, classes[k], classes[i], classes[j])
            end
        end
        push!(M, Hermitian(Mk))
    end
    return ℍVector(M)
end

function dixons_alg(classes, M, h)
    class_orders = []
    for c in classes
        push!(class_orders, length(c))
    end
    g = (2(h)^0.5)
    e = lcm(class_orders...)
    if e == 1
        throw(ArgumentError("I can't handle this yet"))
    end
    primes = get_primes(100) # This will probably need to be adjusted in case the lcm is large
    p = 0
    for prm in primes
        if (prm-1) % e == 0
            if prm > g
                p = prm
                break
            end
        end
    end
    λ = 1:p-1
    println("prime:$p")
    c = 1
    jimbo = nothing
    flag = false
    V = Array{Int64}(undef, length(M), 0)
    println("Doin stuff")
    for k = 1:length(M)
        for l in λ
            pizza = (M[k]-l*I) .% p
            nl = nullspace(pizza)
            sz = size(nl)[2]
            if sz == 0
                continue
            elseif sz == 1
                V = hcat(V, mynorm(nl))
            else
                #jimbo = mapslices(mynorm, nl, dims=(1))
                for v in 1:sz
                    mart = mynorm(nl[:,v])
                    V = hcat(V, mart)
                end
            end
        end
    end
    for v = 1:3
        A = (M[3]-3*I) .% p
        println(A*V[:,v])
    end
end

function mynorm(v)
    b = 0
    for i in v
        if i != 0
            b = i
            break
        end
    end
    return round.(v / b)
end

function get_primes(n)
    primes = []
    for i = 2:n
        if is_prime(i)
            push!(primes, i)
        end
    end
    return primes
end

function is_prime(n)
    if n == 2
        return true
    elseif n == 3
        return true
    end
    sn = floor((n)^0.5)
    for i = 2:sn
        if n % i == 0
            return false
        end
    end
    return true
end