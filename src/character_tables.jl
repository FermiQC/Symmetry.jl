import Base
using LinearAlgebra
using Diagonalizations
include("Shapes.jl")
include("Transformations.jl")
include("PointGroupGenerators.jl")
include("MultiplicationTable.jl")
include("ctab_hard.jl")
tol = 1E-5


function Base.:(==)(A::Symel, B::Symel)
    if sum(A.rrep .- B.rrep) < tol
        return true
    else
        return false
    end
end

function parse_ctfile(ctfilename::String)
    ctstring = read(ctfilename, String)
    parse_ctstring(ctstring)
end

function parse_ctstring(ctstring::String)
    pg_rgx = r"PG:(\w*)"
    ops_rgx = r"Ops:\s*\[([\w,]*)\]"
    irreps_rgx = r"Irreps:\s*\[([\w,]*)\]"
    chars_rgx = r"#Ordered Characters\s*([-\d\s,]*)\s*"
    pg = match(pg_rgx, ctstring).captures[1]
    classes = split(match(ops_rgx, ctstring).captures[1], ",")
    irreps = split(match(irreps_rgx, ctstring).captures[1], ",")
    chars = split(match(chars_rgx, ctstring).captures[1], r"[,\s]+")
    len = size(ops,1)
    char_arr = zeros(len, len)
    for i = 1:len
        for j = 1:len
            char_arr[i,j] = eval(Meta.parse(chars[len*(j-1)+i]))
        end
    end
    symels = ops_to_symels(ops)
    ct = Chartable(pg, symels, irreps, char_arr)
    return ct
end

function pg_to_symels(PG)
    pg = parse_pg_str(PG)
    println(pg)
    symels = [Symel("E", [1 0 0; 0 1 0; 0 0 1])]
    σh = [1 0 0; 0 1 0; 0 0 -1]
    if pg.family == "C"
        if pg.subfamily == "h"
            push!(symels, Symel("sigmah", σh)) # sigmah
            if pg.n % 2 == 0
                push!(symels, Symel("i", i()))
            end
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n)
            symels = vcat(symels, cns, sns)
        elseif pg.subfamily == "v"
            cns = generate_Cn(pg.n)
            if pg.n % 2 == 0
                n = div(pg.n, 2)
                σds = generate_σd(n)
            else
                n = pg.n
                σds = []
            end
            σvs = generate_σv(n)
            symels = vcat(symels, cns, σvs, σds)
        elseif pg.subfamily == "s"
            push!(symels, Symel("sigmah", σh))
        elseif pg.subfamily == "i"
            push!(symels, Symel("i", i()))
        elseif isnothing(pg.subfamily)
            cns = generate_Cn(pg.n)
            symels = vcat(symels, cns)
        else
            throw(ArgumentError("Unidentified Point Group!"))
        end
    elseif pg.family == "D"
        if pg.subfamily == "h"
            push!(symels, Symel("sigmah", σh))
            if pg.n % 2 == 0
                push!(symels, Symel("i", i()))
                n = div(pg.n, 2)
                σds = generate_σd(n)
                c2s = generate_C2(2*pg.n)
            else
                n = pg.n
                σds = []
                c2s = generate_C2(pg.n)
            end
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n)
            σvs = generate_σv(n)
            #c2s = generate_C2(pg.n)
            symels = vcat(symels, cns, sns, σvs, σds, c2s)
        elseif pg.subfamily == "d"
            if pg.n % 2 != 0
                push!(symels, Symel("i", i()))
            end
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n * 2, true)
            σds = generate_σd(pg.n)
            c2s = generate_C2(pg.n)
            symels = vcat(symels, cns, sns, σds, c2s)
        elseif isnothing(pg.subfamily)
            cns = generate_Cn(pg.n)
            c2s = generate_C2(pg.n)
            symels = vcat(symels, cns, c2s)
        else
            throw(ArgumentError("Oh shit, the trout population"))
        end
    elseif pg.family == "S"
        if isnothing(pg.subfamily) & (pg.n % 2 == 0)
            n = div(pg.n, 2)
            if n % 2 != 0
                push!(symels, Symel("i", i()))
            end
            cns = generate_Cn(n)
            sns = generate_Sn(pg.n, true)
            symels = vcat(symels, cns, sns)
        else
            throw(ArgumentError("Oh shit, the trout population"))
        end
    else
        if pg.family == "T"
            if pg.subfamily == "h"
                Ths = generate_Th()
                symels = vcat(symels, Ths)
            elseif pg.subfamily == "d"
                Tds = generate_Td()
                symels = vcat(symels,Tds)
            else
                Ts = generate_T()
                symels = vcat(symels,Ts)
            end
        elseif pg.family == "O"
            if pg.subfamily == "h"
                Ohs = generate_Oh()
                symels = vcat(symels, Ohs)
            else
                Os = generate_O()
                symels = vcat(symels, Os)
            end
        elseif pg.family == "I"
            if pg.subfamily == "h"
                Ihs = generate_Ih()
                symels = vcat(symels, Ihs)
            else
                Is = generate_I()
                symels = vcat(symels, Is)
            end
        else
            throw(ArgumentError("Unidentified Point Group!"))
        end
    end
    return symels
end

function parse_pg_str(s)
    re = r"([A-Z]+)(\d+)?([a-z]+)?"
    m = match(re, s)
    family, n, subfamily = m.captures
    if !isnothing(n)
        n = parse(Int, n)
    end
    if !isnothing(subfamily)
        subfamily = string(subfamily)
    end
    family = string(family)
    return PG(s, family, n, subfamily)
end

function pg_to_chartab(PG)
    pg = parse_pg_str(PG)
    irreps = []
    if pg.family == "C"
        if pg.subfamily == "v"
            println("Cv")
        elseif pg.subfamily == "h"
            println("Ch")
        else
            println("C")
        end
    elseif pg.family == "D"
        if pg.subfamily == "d"
            println("Dd")
        elseif pg.subfamily == "h"
            println("Dh")
        else
            println("D")
        end
    elseif pg.family == "S"
        println("S")
    else
        if pg.family == "T"
            if pg.subfamily == "h"
                println("Th")
            elseif pg.subfamily == "d"
                println("Td")
            else
                println("T")
            end
        elseif pg.family == "O"
            if pg.subfamily == "h"
                println("Oh")
            else
                println("O")
            end
        elseif pg.family == "I"
            if pg.subfamily == "h"
                println("Ih")
            else
                println("I")
            end
        else
            throw(ArgumentError("Unrecognized Point Group"))
        end
    end
end

function char_table_from_str(PG)
    symels = pg_to_symels(PG)
    mtab = build_mult_table(symels)
    println(mtab)
    classes = find_classes(mtab)
    println(classes)
    M = build_Mk(mtab, classes)
    #F = ajd(M; algorithm=:JADE)
    #println(F.ev)
    #G = zeros(length(classes), length(classes))
    #for i = 1:length(classes)
    #    G[:,i] += normalize(F.F[:,i])
    #end
    #println(G)
    ##dixons_alg(classes, M, length(symels))
end