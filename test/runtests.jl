using DiscreteMarkovChains
using Test
using SafeTestsets

@time begin
    @time @safetestset "Utils" begin include("Utils.jl") end
    @time @safetestset "MCCore" begin include("MCCore.jl") end
end
