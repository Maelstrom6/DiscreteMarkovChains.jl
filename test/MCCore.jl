using DiscreteMarkovChains
using Test

@testset "Expected Failures" begin
    # Different lengths
    @test_throws(
        ErrorException(
            "The state space, [1, 2], and "*
            "transition matrix should be the same size."
        ),
        DiscreteMarkovChain([1, 2], [1]),
    )

    # Repeated elements in the state space
    @test_throws(
        ErrorException("The state space, [1, 1], must have unique elements."),
        DiscreteMarkovChain([1, 1], [1 0; 0 1]),
    )

    # Not row stochastic
    @test_throws(
        ErrorException(
            "The transition matrix, [0.5 0.6; 1.0 0.0], "*
            "should be row-stochastic (each row must sum up to 1)."
        ),
        DiscreteMarkovChain([1, 2], [0.5 0.6; 1.0 0.0]),
    )
end

@testset "Communication Classes" begin
    # Empty test
    T = [][:,:]
    X = DiscreteMarkovChain([], T)
    @test periodicities(X) == ([], [], [])

    state_space = [1, 2, 3]
    T = [
        0 1 0;
        1 0 0;
        0 0 1;
    ]
    X = DiscreteMarkovChain(state_space, T)
    @test periodicities(X) == ([[1, 2], [3]], [true, true], [2, 1])

    T = [
        0 5 5 0 0;
        0 0 0 10 0;
        5 0 5 0 0;
        0 10 0 0 0;
        0 3 0 3 4;
    ]/10
    X = DiscreteMarkovChain(0:4, T)
    classes, recurrence, periods = periodicities(X)
    @test classes == [[1, 3], [0, 2], [4]]
    @test recurrence == [true, false, false]
    @test periods == [2, 1, 1]

    T = [
        0 0 0 10 0 0;
        5 0 5 0 0 0;
        0 4 0 0 0 6;
        10 0 0 0 0 0;
        0 10 0 0 0 0;
        0 0 0 5 5 0;
    ]/10
    X = DiscreteMarkovChain(0:5, T)
    classes, recurrence, periods = periodicities(X)
    @test classes == [[0, 3], [1, 2, 5, 4]]
    @test recurrence == [true, false]
    @test periods == [2, 2]

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
    ]/10
    X = DiscreteMarkovChain(0:9, T)
    classes, recurrence, periods = periodicities(X)
    @test classes == [[0, 3, 6, 7], [1], [2, 8, 9], [5], [4]]
    @test recurrence == [true, true, false, true, false]
    @test periods == [1, 1, 1, 1, 1]
end

@testset "Decompose and Canonical Form" begin
    # Empty test
    T = Array{Any}(undef, 0, 0)
    X = DiscreteMarkovChain([], T)
    @test decompose(X) == ([], T, T, T)

    T = [
        1 0 0;
        1//3 1//3 1//3;
        0 1//4 3//4;
    ]
    X = DiscreteMarkovChain(["Sunny", "Cloudy", "Rainy"], T)
    @test decompose(X) == (
        ["Sunny", "Cloudy", "Rainy"],
        T[1:1, 1:1], T[2:3, 1:1], T[2:3, 2:3],
    )
    @test canonical_form(X) == (["Sunny", "Cloudy", "Rainy"], T)

    T = [
        1//4 3//4 0;
        1//3 1//3 1//3;
        0 1//4 3//4;
    ]
    X = DiscreteMarkovChain([1, 2, 3], T)
    @test decompose(X) == ([1, 2, 3], T, Array{Any}(undef, 0, 3), Array{Any}(undef, 0, 0))
    @test canonical_form(X) == ([1, 2, 3], T)

    T = [
        1 0 0 0 0;
        1//2 0 1//2 0 0;
        0 1//2 0 1//2 0;
        0 0 1//2 0 1//2;
        0 0 0 0 1;
    ]
    X = DiscreteMarkovChain(0:4, T)
    states, A, B, C = decompose(X)
    @test states == [0, 4, 1, 2, 3]
    @test A == [1 0; 0 1]
    @test B == [1//2 0; 0 0; 0 1//2]
    @test C == [0 1//2 0; 1//2 0 1//2; 0 1//2 0]
    states, new_matrix = canonical_form(X)
    @test states == [0, 4, 1, 2, 3]
    desired_new_matrix = [
        1 0 0 0 0;
        0 1 0 0 0;
        1//2 0 0 1//2 0;
        0 0 1//2 0 1//2;
        0 1//2 0 1//2 0;
    ]
    @test new_matrix == desired_new_matrix
end

@testset "Regular, Ergodic, Absorbing" begin
    # Empty test
    T = Array{Any}(undef, 0, 0)
    X = DiscreteMarkovChain([], T)
    @test !is_regular(X)
    @test !is_ergodic(X)
    @test !is_absorbing(X)

    T = [
        0 4 0 0 0;
        1 0 3 0 0;
        0 2 0 2 0;
        0 0 3 0 1;
        0 0 0 4 0;
    ]/4.0
    X = DiscreteMarkovChain(1:5, T)
    @test !is_regular(X)
    @test is_ergodic(X)
    @test !is_absorbing(X)

    T = [
        0 1;
        1 0;
    ]
    X = DiscreteMarkovChain(1:2, T)
    @test !is_regular(X)
    @test is_ergodic(X)
    @test !is_absorbing(X)

    T = [
        2 1 1;
        2 0 2;
        1 1 2;
    ]/4
    X = DiscreteMarkovChain(1:3, T)
    @test is_regular(X)
    @test is_ergodic(X)
    @test !is_absorbing(X)

    T = [
        1 1;
        1 1;
    ]/2
    X = DiscreteMarkovChain(1:2, T)
    @test is_regular(X)
    @test is_ergodic(X)
    @test !is_absorbing(X)
end

@testset "Stationary Distribution" begin
    # Empty test
    T = Array{Any}(undef, 0, 0)
    X = DiscreteMarkovChain([], T)
    @test stationary_distribution(X) == Any[]

    T = [
        1//5 2//5 2//5;
        1//10 1//2 2//5;
        3//5 3//10 1//10;
    ]
    X = DiscreteMarkovChain(1:3, T)
    @test stationary_distribution(X) == [11//39, 16//39, 4//13]

    T = [
        0 1 0 0 0;
        .25 0 .75 0 0;
        0 .5 0 .5 0;
        0 0 .75 0 .25;
        0 0 0 1 0;
    ]
    X = DiscreteMarkovChain(0:4, T)
    @test stationary_distribution(X) ≈ [.0625, .2500, .3750, .2500, .0625]
end

@testset "Fundamental Matrix" begin
    # Empty test
    T = Array{Any}(undef, 0, 0)
    X = DiscreteMarkovChain([], T)
    @test fundamental_matrix(X) == Array{Any}(undef, 0, 0)

    T = [
        1 0 0 0 0;
        1//2 0 1//2 0 0;
        0 1//2 0 1//2 0;
        0 0 1//2 0 1//2;
        0 0 0 0 1;
    ]
    X = DiscreteMarkovChain([0, 1, 2, 3, 4], T)
    desired_new_matrix = [
        3//2 1 1//2;
        1 2 1;
        1//2 1 3//2;
    ]
    @test fundamental_matrix(X) == desired_new_matrix
end

@testset "Expected Time To Absorption" begin
    # Empty test
    T = Array{Any}(undef, 0, 0)
    X = DiscreteMarkovChain([], T)
    @test mean_time_to_absorption(X) == Any[]

    T = [
        1 0 0 0 0;
        1//2 0 1//2 0 0;
        0 1//2 0 1//2 0;
        0 0 1//2 0 1//2;
        0 0 0 0 1;
    ]
    X = DiscreteMarkovChain([0, 1, 2, 3, 4], T)
    @test mean_time_to_absorption(X) == [2, 3, 2]

    T = [
        1 0 0 0 0;
        0 1 0 0 0;
        0 0 0 1 0;
        1//4 1//4 1//4 0 1//4;
        1//2 0 0 1//2 0;
    ]
    X = DiscreteMarkovChain([0, 1, 2, 3, 4], T)
    @test mean_time_to_absorption(X) == [12//5, 7//5, 6//5]
end

@testset "Exit Probability" begin
    # Empty test
    T = Array{Any}(undef, 0, 0)
    X = DiscreteMarkovChain([], T)
    @test exit_probabilities(X) == Array{Any}(undef, 0, 0)

    T = [
        1 0 0 0 0;
        0 1 0 0 0;
        0 0 0 1 0;
        1//4 1//4 1//4 0 1//4;
        1//2 0 0 1//2 0;
    ]
    X = DiscreteMarkovChain([0, 1, 2, 3, 4], T)
    desired_new_matrix = [
        3//5 2//5;
        3//5 2//5;
        4//5 1//5;
    ]
    @test exit_probabilities(X) == desired_new_matrix

    T = [
        2//5 3//5 0 0 0 0;
        9//10 1//10 0 0 0 0;
        0 0 1 0 0 0;
        0 0 0 0 1 0;
        0 1//4 1//4 1//4 0 1//4;
        0 1//2 0 0 1//2 0;
    ]
    X = DiscreteMarkovChain([0, 1, 2, 3, 4, 5], T)
    desired_new_matrix = [
        0 6//10 4//10;
        0 6//10 4//10;
        0 8//10 2//10;
    ]
    @test exit_probabilities(X) == desired_new_matrix
end

@testset "First Passage" begin
    # Empty test
    T = Array{Any}(undef, 0, 0)
    X = DiscreteMarkovChain([], T)
    @test first_passage_probabilities(X, 1) == Array{Any}(undef, 0, 0)

    T = [
        2//10 4//10 4//10;
        3//10 3//10 4//10;
        5//10 4//10 1//10;
    ]
    X = DiscreteMarkovChain(["A", "B", "C"], T)
    @test first_passage_probabilities(X, 1) == T

    @test first_passage_probabilities(X, 1, "A", "C") == 4//10
    @test first_passage_probabilities(X, 2, "A", "C") == 24//100
    @test first_passage_probabilities(X, 3, "A", "C") == 144//1000
    @test first_passage_probabilities(X, 4, "A", "C") == 864//10000

    @test first_passage_probabilities(X, 1, "A", "A") == 2//10
    @test first_passage_probabilities(X, 2, "A", "A") == 32//100
    @test first_passage_probabilities(X, 3, "A", "A") == 184//1000
    @test first_passage_probabilities(X, 4, "A", "A") == 1152//10000
end

@testset "Mean Recurrence Time" begin
    # Empty test
    T = Array{Any}(undef, 0, 0)
    X = DiscreteMarkovChain([], T)
    @test mean_recurrence_time(X) == Any[]

    T = [
        0 1 0;
        0 0 1;
        1 0 0;
    ]
    X = DiscreteMarkovChain(T)
    @test mean_recurrence_time(X) ≈ [3, 3, 3]

    T = [
        0 1 0 0 0;
        .25 0 .75 0 0;
        0 .5 0 .5 0;
        0 0 .75 0 .25;
        0 0 0 1 0;
    ]
    X = DiscreteMarkovChain(0:4, T)
    @test mean_recurrence_time(X) ≈ [16.0000, 4.0000, 8/3, 4.0000, 16.000]
end

@testset "Mean First Passage Time" begin
    # Empty test
    T = Array{Any}(undef, 0, 0)
    X = DiscreteMarkovChain([], T)
    @test mean_first_passage_time(X) == Array{Any}(undef, 0, 0)

    T = [
        1//2 1//4 1//4;
        1//2 0 1//2;
        1//4 1//4 1//2;
    ]
    X = DiscreteMarkovChain(["A", "B", "C"], T)
    desired_new_matrix = [
        0 4 10//3;
        8//3 0 8//3;
        10//3 4 0;
    ]
    @test mean_first_passage_time(X) == desired_new_matrix

    T = [
        0 1 0 0 0;
        .25 0 .75 0 0;
        0 .5 0 .5 0;
        0 0 .75 0 .25;
        0 0 0 1 0;
    ]
    X = DiscreteMarkovChain(0:4, T)
    desired_new_matrix = [
        0 1 8/3 19/3 64/3;
        15 0 5/3 16/3 61/3;
        56/3 11/3 0 11/3 56/3;
        61/3 16/3 5/3 0 15;
        64/3 19/3 8/3 1 0;
    ]
    @test mean_first_passage_time(X) ≈ desired_new_matrix
end
