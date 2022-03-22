@testset "Symmetry Elements" begin
    for t in allofem
        @test Symmetry.pg_to_symels(t) == eval(Meta.parse(t*"s"))
    end
end