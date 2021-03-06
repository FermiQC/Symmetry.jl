module Symmetry
using LinearAlgebra
import Base

tol = 1E-5

include("Shapes.jl")
include("Transformations.jl")
include("PointGroupGenerators.jl")
include("ctab_hard.jl")
include("character_tables.jl")

end