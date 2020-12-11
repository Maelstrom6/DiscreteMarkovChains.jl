abstract type AbstractDiscreteMarkovChain <: AbstractMarkovChain end

"""
    DiscreteMarkovChain(transition_matrix)
    DiscreteMarkovChain(state_space, transition_matrix)
    DiscreteMarkovChain(continuous_markov_chain)

Creates a new discrete Markov chain object.

# Arguments
- `state_space`: The names of the states that make up the Markov chain.
- `transition_matrix`: The single step transition probability matrix.
- `continuous_markov_chain`: An instance of `ContinuousMarkovChain`.

# Examples
The following shows a basic Sunny-Cloudy-Rainy weather model.
```jldoctest DiscreteMarkovChain
using DiscreteMarkovChains
T = [
    0.9 0.1 0;
    0.5 0.2 0.3;
    0.1 0.4 0.5
]
X = DiscreteMarkovChain(["Sunny", "Cloudy", "Rainy"], T)
println(state_space(X))

# output

["Sunny", "Cloudy", "Rainy"]
```

```jldoctest DiscreteMarkovChain
println(transition_matrix(X))

# output

[0.9 0.1 0.0; 0.5 0.2 0.3; 0.1 0.4 0.5]
```

# References
1. [Wikipedia](https://en.wikipedia.org/wiki/Markov_chain#Discrete-time_Markov_chain)
2. [Dartmouth College](https://www.dartmouth.edu/~chance/teaching_aids/books_articles/probability_book/Chapter11.pdf)
"""
struct DiscreteMarkovChain <: AbstractDiscreteMarkovChain
    state_space
    transition_matrix
    function DiscreteMarkovChain(state_space, transition_matrix)
        check(state_space, transition_matrix, DiscreteMarkovChain)
        new(state_space, transition_matrix)
    end
end
"""
    required_row_sum(type)

# Arguments
- `type`: The type of Markov chain. It can be
`DiscreteMarkovChain` or `ContinuousMarkovChain`.

# Returns
The number that each row in the transition matrix should sum up to.
"""
required_row_sum(::Core.Type{<:AbstractDiscreteMarkovChain}) = 1
function DiscreteMarkovChain(transition_matrix)
    return DiscreteMarkovChain(1:(size(transition_matrix)[1]), transition_matrix)
end

"""
    characteristic_matrix(::AbstractDiscreteMarkovChain)
    characteristic_matrix(::AbstractContinuousMarkovChain)

# Definitions
Many derivations and interesting ideas about Markov chains
involve the identity matrix or zero matrix somewhere along
the line. Most of the time, the identity matrix appears
more often in discrete Markov chains. Instead of the identity
matrix, the zero matrix appears in its place for continuous
Markov chains.

# Returns
The identity matrix if its argument is an instance of
`AbstractDiscreteMarkovChain`. The zero matrix if its
argument is an instance of `AbstractContinuousMarkovChain`
"""
characteristic_matrix(::AbstractDiscreteMarkovChain) = LinearAlgebra.I

"""
    periodicities(x::AbstractDiscreteMarkovChain)

A more advanced version of `communication_classes`
designed for discrete Markov chains. It is the same as
`communication_classes` but it returns periodicities as well.

# Definitions
The period, ``d_i`` of a state ``i`` is the greatest common denominator
of all integers ``n ∈ \\mathbb{N}`` for which ``p^{(n)}_{i,i} > 0``.
Written more succinctly,
```math
d_i = \\text{gcd}\\{ n ∈ \\mathbb{N} | p^{(n)}_{i,i} > 0 \\}
```
If ``d_i=1`` then state ``i`` is said to be aperiodic.

# Arguments
- `x`: some kind of discrete Markov chain.

# Returns
A tuple containing 3 arrays.
- This first array contains C arrays which store the states that communicate.
- The second array is an array of Bool where the ith value is true if the
  ith communication class is recurrent.
- The third array is the periodicity of each communication class.

# Examples
```jldoctest periodicities
using DiscreteMarkovChains
T = [
    1 0;
    0 1;
]
X = DiscreteMarkovChain(T)

periodicities(X)

# output

([[1], [2]], Any[true, true], Any[1, 1])
```

So the Markov chain has two communication classes and both are recurrent and aperiodic.

```jldoctest periodicities
T = [
    0.0 1.0 0.0;
    1.0 0.0 0.0;
    0.1 0.2 0.7;
]
X = DiscreteMarkovChain(["Sunny", "Cloudy", "Rainy"], T)

periodicities(X)

# output

([["Sunny", "Cloudy"], ["Rainy"]], Any[true, false], Any[2, 1])
```

So the Sunny and Cloudy states communicate and are recurrent with period 2.
The Rainy state is transient and aperiodic.
Note that this is not a very good weather model since once it
stops raining, it will never rain again. Also, each day after that,
the process will oscillate between Sunny and Cloudy.
"""
function periodicities(x::AbstractDiscreteMarkovChain)
    SI = state_index(x)
    T = transition_matrix(x)

    classes, recurrence = communication_classes(x)
    index_classes = [[SI[state] for state in class] for class in classes]

    periods = []
    for class in index_classes
        submatrix = T[class, class]
        push!(periods, breadth_first_search(submatrix))
    end
    return classes, recurrence, periods
end

"""
    is_regular(x)

# Definitions
A Markov chain is called a regular chain if some power of the
transition matrix has only positive elements. This is equivalent
to being ergodic and aperiodic.

# Arguments
- `x`: some kind of discrete Markov chain.

# Returns
`true` if the Markov chain, `x`, is regular.

# Examples
We will set up a matrix with 2 communication classes and
show that it is not regular.

```jldoctest is_regular
using DiscreteMarkovChains
T = [
    0 1 0;
    1 0 0;
    0 0 1;
]
X = DiscreteMarkovChain(T)

is_regular(X)

# output

false
```

Repeat the above but now all states communicate.

```jldoctest is_regular
T = [
    0.0 0.5 0.5;
    0.0 0.0 1.0;
    1.0 0.0 0.0;
]
X = DiscreteMarkovChain(T)

is_regular(X)

# output

true
```

Notice how a periodic chain is not regular even though
there is only one communication class.

```jldoctest is_regular
T = [
    0 1 0;
    0 0 1;
    1 0 0;
]
X = DiscreteMarkovChain(T)

is_regular(X)

# output

false
```
"""
function is_regular(x::AbstractDiscreteMarkovChain)
    classes, _, periods = periodicities(x)
    if length(classes) == 0
        return false
    end
    return (length(classes) == 1) & (periods[1] == 1)
end

"""
    first_passage_probabilities(x, t, i=missing, j=missing)

# Definitions
This is the probability that the process enters state ``j``
for the first time at time ``t`` given that the process started
in state ``i`` at time 0. That is, ``f^{(t)}_{i,j}``. If no `i`
or `j` is given, then it will return a matrix instead with
entries ``f^{(t)}_{i,j}`` for `i` and `j` in the state space of `x`.

# Why Do We Use A Slow Algorithm?
So that `t` can be symbolic if nessesary. That is, if symbolic math
libraries want to use this library, it will pose no hassle.

# Arguments
- `x`: some kind of Markov chain.
- `t`: the time to calculate the first passage probability.
- `i`: the state that the prcess starts in.
- `j`: the state that the process must reach for the first time.

# Returns
A scalar value or a matrix depending on whether `i` and `j` are given.

# Examples
```jldoctest first_passage_probabilities
using DiscreteMarkovChains
T = [
    0.1 0.9;
    0.3 0.7;
]
X = DiscreteMarkovChain(T)

first_passage_probabilities(X, 2)

# output

2×2 Array{Float64,2}:
 0.27  0.09
 0.21  0.27
```

If `X` has a custom state space,
then `i` and `j` must be in that state space.

```jldoctest first_passage_probabilities
T = [
    0.1 0.9;
    0.3 0.7;
]
X = DiscreteMarkovChain(["Sunny", "Rainy"], T)

first_passage_probabilities(X, 2, "Sunny", "Rainy")

# output

0.09000000000000001
```

Notice how this is the (1, 2) entry in the first example.

# References
1. [University of Windsor](https://scholar.uwindsor.ca/cgi/viewcontent.cgi?article=1125&context=major-papers)
2. [Durham University](http://maths.dur.ac.uk/stats/courses/ProbMC2H/_files/handouts/1516MarkovChains2H.pdf)
"""
function first_passage_probabilities(
    x::AbstractDiscreteMarkovChain, t, i=missing, j=missing
)
    S = state_space(x)
    T = transition_matrix(x)
    n = length(S)

    if n == 0
        return transition_matrix(x)
    end

    js = 1:n  # The columns to loop through
    calc_i_ne_j = true  # Calculate the off-diagonals
    calc_i_eq_j = true  # Calculate the diagonals
    if (i !== missing) && (j !== missing)
        i = state_index(x)[i]
        j = state_index(x)[j]
        js = [j]
        if i == j
            calc_i_ne_j = false
        else
            calc_i_eq_j = false
        end
    end

    Ft = zeros(eltype(T), n, n)  # Empty matrix

    # If i != j
    if calc_i_ne_j
        for j in js

            P0 = copy(T)
            P0[1:n, j] = zeros(eltype(T), n, 1)
            F = P0^(t-1) * T
            Ft[1:n, j] = F[1:n, j]
        end
    end

    # If i == j
    if calc_i_eq_j
        for j in js

            P_ = copy(T)
            P_[j, 1:n] = zeros(eltype(T), 1, n)

            Pnew = zeros(eltype(T), 2*n, 2*n)
            Pnew[1:n, 1:n] = T
            Pnew[(n+1):(2*n), (n+1):(2*n)] = P_
            Pnew[n+j, 1:n] = T[j, 1:n]

            P0 = copy(Pnew)
            P0[1:(2*n), j] = zeros(eltype(T), 2*n, 1)

            F = P0^(t - 1) * Pnew

            Ft[j, j] = F[n+j, j]
        end
    end

    if (i !== missing) && (j !== missing)
        return Ft[i, j]
    end
    return Ft
end
