using GasphaseReactions
using Test

@testset "GasphaseReactions.jl" begin
    retcode = GasphaseReactions.run("gaschem/gaschem.xml","lib/")
    @test abs(retcode) < 1e-9
end
