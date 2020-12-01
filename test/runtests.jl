using DiscreteMarkovChains
using Test
using SafeTestsets

@time begin
    @time @safetestset "Utils" begin include("Utils.jl") end
    @time @safetestset "Core" begin include("Core.jl") end
end
