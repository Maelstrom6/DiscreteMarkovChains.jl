using DiscreteMarkovChains
using Test
using LinearAlgebra

@testset "Converting To Discrete" begin
    # Empty test
    T = [][:,:]
    X = ContinuousMarkovChain(T)
    @test DiscreteMarkovChain(X) == DiscreteMarkovChain(T)

    state_space = [1, 2, 3]
    T = [
        0 1 0;
        1 0 0;
        0 0 1;
    ] - I
    X = ContinuousMarkovChain(state_space, T)

    T1 = transition_matrix(DiscreteMarkovChain(X))
    T2 = transition_matrix(DiscreteMarkovChain(state_space, exp(T)))
    @test T1 == T2
end

@testset "Converting To Continuous" begin
    # Empty test
    T = [][:,:]
    X = DiscreteMarkovChain(T)
    @test ContinuousMarkovChain(X) == ContinuousMarkovChain(T)

    state_space = [1, 2, 3]
    T = [
        0.2 0.8 0;
        0 0.2 0.8;
        0 0 1;
    ]
    X = DiscreteMarkovChain(state_space, T)

    Q1 = transition_matrix(ContinuousMarkovChain(X))
    Q2 = transition_matrix(ContinuousMarkovChain(state_space, log(T)))
    @test Q1 == Q2
end

@testset "Embedded Chains" begin
    # Empty test
    T = [][:,:]
    X = ContinuousMarkovChain(T)
    @test embedded(X) == DiscreteMarkovChain(T)

    X = DiscreteMarkovChain(T)
    @test embedded(X) == X

    T = [
        -2 1 1;
        1 -1 0;
        2 3 -5;
    ]
    X = ContinuousMarkovChain(T)
    desired_new_matrix = [
        0 0.5 0.5;
        1 0 0;
        0.4 0.6 0;
    ]
    @test transition_matrix(embedded(X)) ≈ desired_new_matrix
end

@testset "Communication Classes" begin
    # Empty test
    T = [][:,:]
    X = ContinuousMarkovChain(T)
    @test communication_classes(X) == ([], [])

    state_space = [1, 2, 3]
    T = [
        0 1 0;
        1 0 0;
        0 0 1;
    ] - I
    X = ContinuousMarkovChain(state_space, T)
    @test communication_classes(X) == ([[1, 2], [3]], [true, true])

    T = [
        0 5 5 0 0;
        0 0 0 10 0;
        5 0 5 0 0;
        0 10 0 0 0;
        0 3 0 3 4;
    ]/10 - I
    X = ContinuousMarkovChain(0:4, T)
    classes, recurrence = communication_classes(X)
    @test classes == [[1, 3], [0, 2], [4]]
    @test recurrence == [true, false, false]

    T = [
        0 0 0 10 0 0;
        5 0 5 0 0 0;
        0 4 0 0 0 6;
        10 0 0 0 0 0;
        0 10 0 0 0 0;
        0 0 0 5 5 0;
    ]/10 - I
    X = ContinuousMarkovChain(0:5, T)
    classes, recurrence = communication_classes(X)
    @test classes == [[0, 3], [1, 2, 5, 4]]
    @test recurrence == [true, false]

    T = [
        2 0 0 3 0 0 3 2 0 0;
        0 10 0 0 0 0 0 0 0 0;
        0 2 2 0 0 0 0 0 3 3;
        0 0 0 3 0 0 6 1 0 0;
        0 0 0 0 5 5 0 0 0 0;
        0 0 0 0 0 10 0 0 0 0;
        4 0 0 5 0 0 1 0 0 0;
        2 0 0 4 0 0 2 2 0 0;
        3 0 1 0 0 0 0 0 4 2;
        0 0 4 0 0 0 0 0 3 3;
    ]/10 - I
    X = ContinuousMarkovChain(0:9, T)
    classes, recurrence = communication_classes(X)
    @test classes == [[0, 3, 6, 7], [1], [2, 8, 9], [5], [4]]
    @test recurrence == [true, true, false, true, false]
end

@testset "Decompose and Canonical Form" begin
    # Empty test
    T = Array{Any}(undef, 0, 0)
    X = ContinuousMarkovChain([], T)
    @test decompose(X) == ([], T, T, T)

    T = [
        1 0 0;
        1//3 1//3 1//3;
        0 1//4 3//4;
    ] - I
    X = ContinuousMarkovChain(["Sunny", "Cloudy", "Rainy"], T)
    @test decompose(X) == (
        ["Sunny", "Cloudy", "Rainy"],
        T[1:1, 1:1], T[2:3, 1:1], T[2:3, 2:3],
    )
    @test canonical_form(X) == (["Sunny", "Cloudy", "Rainy"], T)

    T = [
        1//4 3//4 0;
        1//3 1//3 1//3;
        0 1//4 3//4;
    ] - I
    X = ContinuousMarkovChain([1, 2, 3], T)
    @test decompose(X) == ([1, 2, 3], T, Array{Any}(undef, 0, 3), Array{Any}(undef, 0, 0))
    @test canonical_form(X) == ([1, 2, 3], T)

    T = [
        1 0 0 0 0;
        1//2 0 1//2 0 0;
        0 1//2 0 1//2 0;
        0 0 1//2 0 1//2;
        0 0 0 0 1;
    ] - I
    X = ContinuousMarkovChain(0:4, T)
    states, A, B, C = decompose(X)
    @test states == [0, 4, 1, 2, 3]
    @test A == [1 0; 0 1] - I
    @test B == [1//2 0; 0 0; 0 1//2]
    @test C == [0 1//2 0; 1//2 0 1//2; 0 1//2 0] - I
    states, new_matrix = canonical_form(X)
    @test states == [0, 4, 1, 2, 3]
    desired_new_matrix = [
        1 0 0 0 0;
        0 1 0 0 0;
        1//2 0 0 1//2 0;
        0 0 1//2 0 1//2;
        0 1//2 0 1//2 0;
    ] - I
    @test new_matrix == desired_new_matrix
end

@testset "Regular, Ergodic, Absorbing" begin
    # Empty test
    T = Array{Any}(undef, 0, 0)
    X = ContinuousMarkovChain([], T)
    @test !is_ergodic(X)
    @test !is_absorbing(X)

    T = [
        0 4 0 0 0;
        1 0 3 0 0;
        0 2 0 2 0;
        0 0 3 0 1;
        0 0 0 4 0;
    ]/4.0 - I
    X = ContinuousMarkovChain(1:5, T)
    @test is_ergodic(X)
    @test !is_absorbing(X)

    T = [
        0 1;
        1 0;
    ] - I
    X = ContinuousMarkovChain(1:2, T)
    @test is_ergodic(X)
    @test !is_absorbing(X)

    T = [
        2 1 1;
        2 0 2;
        1 1 2;
    ]/4 - I
    X = ContinuousMarkovChain(1:3, T)
    @test is_ergodic(X)
    @test !is_absorbing(X)

    T = [
        1 1;
        1 1;
    ]/2 - I
    X = ContinuousMarkovChain(1:2, T)
    @test is_ergodic(X)
    @test !is_absorbing(X)
end

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

@testset "Fundamental Matrix" begin
    # Empty test
    T = Array{Any}(undef, 0, 0)
    X = ContinuousMarkovChain([], T)
    @test fundamental_matrix(X) == Array{Any}(undef, 0, 0)

    # http://www.mat.ufrgs.br/~giacomo/Softwares/R/Dobrow%20-%20Stochastic%20with%20R/ch7.pdf
    # Example 7.16
    T = [
        0 0 0;
        1 -1 0;
        1 1 -2;
    ]
    X = ContinuousMarkovChain(T)
    desired_new_matrix = [
        1 0;
        1/2 1/2;
    ]
    @test fundamental_matrix(X) == desired_new_matrix
end

@testset "Expected Time To Absorption" begin
    # Empty test
    T = Array{Any}(undef, 0, 0)
    X = ContinuousMarkovChain([], T)
    @test mean_time_to_absorption(X) == Any[]

    T = [
        0 0 0;
        4 -4 0;
        0 2 -2;
    ]
    X = ContinuousMarkovChain(T)
    @test mean_time_to_absorption(X) == [0.25, 0.75]
end

@testset "Exit Probability" begin
    # Empty test
    T = Array{Any}(undef, 0, 0)
    X = ContinuousMarkovChain([], T)
    @test exit_probabilities(X) == Array{Any}(undef, 0, 0)

    T = [
        0 0 0;
        0 -1 1;
        2 2 -4;
    ]
    X = ContinuousMarkovChain(T)
    @test exit_probabilities(X) == [1; 1][:, :]

    T = [
        0 0 0;
        0 0 0;
        8 8 -16;
    ]
    X = ContinuousMarkovChain(T)
    @test exit_probabilities(X) == [0.5 0.5]
end

@testset "Mean Recurrence Time" begin
    # Empty test
    T = Array{Any}(undef, 0, 0)
    X = ContinuousMarkovChain([], T)
    @test mean_recurrence_time(X) == Any[]

    T = [
        -1 1 0;
        0 -1 1;
        1 0 -1;
    ]
    X = ContinuousMarkovChain(T)
    @test mean_recurrence_time(X) == [3, 3, 3]

    T = [
        -1 1 0;
        0 -1 1;
        2 2 -4;
    ]
    X = ContinuousMarkovChain(T)
    @test mean_recurrence_time(X) == [3.5, 1.75, 1.75]
end

@testset "Mean First Passage Time" begin
    # Empty test
    T = Array{Any}(undef, 0, 0)
    X = ContinuousMarkovChain([], T)
    @test mean_first_passage_time(X) == Array{Any}(undef, 0, 0)

    T = [
        -1 1 0;
        0 -1 1;
        1 0 -1;
    ]
    X = ContinuousMarkovChain(T)
    @test mean_first_passage_time(X) == [
        0.0 1.0 2.0;
        2.0 0.0 1.0;
        1.0 2.0 0.0;
    ]

    T = [
        -1 1 0;
        0 -1 1;
        2 2 -4;
    ]
    X = ContinuousMarkovChain(T)
    @test mean_first_passage_time(X) == [
        0.0 1.0 2.0;
        2.5 0.0 1.0;
        1.5 0.75 0.0;
    ]
end
