using DiscreteMarkovChains
using Test

@testset "Stationary distribution" begin
    # Empty test
    T = Array{Any}(undef, 0, 0)
    X = ContinuousMarkovChain([], T)
    @test stationary_distribution(X) == Any[]

    T = [-2 2 0; 0 -4 4; 0 0 0]
    X = ContinuousMarkovChain(T)
    @test stationary_distribution(X) == [0, 0, 1]

    T = [
        -25 20 5;
        300 -500 200;
        20 400 -420;
    ]//1000
    X = ContinuousMarkovChain(T)
    @test stationary_distribution(X) == [100//113, 8//113, 5//113]

end
