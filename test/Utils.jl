import DiscreteMarkovChains
using Test

@testset "Strongly Connected Components" begin
    V = [1, 2, 3]
    E = [(1, 2), (2, 1), (3, 1)]
    components = DiscreteMarkovChains.strongly_connected_components(V, E)
    @test components == [[1, 2], [3]]
end

@testset "Breadth First Search" begin
    T = [1][:,:]
    period = DiscreteMarkovChains.breadth_first_search(T)
    @test period == 1

    T = [
        0 1;
        1 0
    ]
    period = DiscreteMarkovChains.breadth_first_search(T)
    @test period == 2

    T = [
        0 0.5 0.5;
        0.5 0.5 0;
        0.5 0.5 0;
    ]
    period = DiscreteMarkovChains.breadth_first_search(T)
    @test period == 1

    # this was a crazy bug
    T = [
        0 4 0 0 0;
        1 0 3 0 0;
        0 2 0 2 0;
        0 0 3 0 1;
        0 0 0 4 0;
    ]
    period = DiscreteMarkovChains.breadth_first_search(T)
    @test period == 2
end
