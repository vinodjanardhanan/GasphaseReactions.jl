using GasphaseReactions
using Test

@testset "GasphaseReactions.jl" begin
    if Sys.isapple() || Sys.islinux()
        retcode = GasphaseReactions.run("gaschem/gaschem.xml","lib/")
        @test abs(retcode) < 1e-9
    elseif Sys.iswindows() 
        retcode = GasphaseReactions.run("gaschem\\gaschem.xml","lib\\")
        @test abs(retcode) < 1e-9
    end
end
