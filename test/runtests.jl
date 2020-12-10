using DiscreteMarkovChains
using Test
using SafeTestsets

@time begin
    @time @safetestset "Utils" begin include("Utils.jl") end
    @time @safetestset "DMC" begin include("MCCore.jl") end
    @time @safetestset "CMC" begin include("CMC.jl") end
end
