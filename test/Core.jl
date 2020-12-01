using DiscreteMarkovChains
using Test


@testset "Communication Classes" begin
    # empty test
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
    # empty test
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
    # empty test
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
    # empty test
    T = Array{Any}(undef, 0, 0)
    X = DiscreteMarkovChain([], T)
    @test stationary_distribution(X) == Array{Any}(undef, 0, 0)

    T = [
        1//5 2//5 2//5;
        1//10 1//2 2//5;
        3//5 3//10 1//10;
    ]
    X = DiscreteMarkovChain(1:3, T)
    @test stationary_distribution(X) == [11//39, 16//39, 4//13]
end

@testset "Fundamental Matrix" begin
    # empty test
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
    # empty test
    T = Array{Any}(undef, 0, 0)
    X = DiscreteMarkovChain([], T)
    @test expected_time_to_absorption(X) == Any[]

    T = [
        1 0 0 0 0;
        1//2 0 1//2 0 0;
        0 1//2 0 1//2 0;
        0 0 1//2 0 1//2;
        0 0 0 0 1;
    ]
    X = DiscreteMarkovChain([0, 1, 2, 3, 4], T)
    @test expected_time_to_absorption(X) == [2, 3, 2]

    T = [
        1 0 0 0 0;
        0 1 0 0 0;
        0 0 0 1 0;
        1//4 1//4 1//4 0 1//4;
        1//2 0 0 1//2 0;
    ]
    X = DiscreteMarkovChain([0, 1, 2, 3, 4], T)
    @test expected_time_to_absorption(X) == [12//5, 7//5, 6//5]
end

@testset "Exit Probability" begin
    # empty test
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
    # empty test
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
