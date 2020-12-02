module DiscreteMarkovChains
import LinearAlgebra

include("Utils.jl")
export DiscreteMarkovChain, state_space, transition_matrix, 
communication_classes, periodicities, decompose, canonical_form,
is_regular, is_ergodic, is_absorbing,
stationary_distribution, fundamental_matrix, 
expected_time_to_absorption, exit_probabilities, first_passage_probabilities,
mean_recurrence_time, mean_first_passage_time

abstract type AbstractMarkovChain end
abstract type AbstractDiscreteMarkovChain <: AbstractMarkovChain end

"""
    DiscreteMarkovChain(transition_matrix)
    DiscreteMarkovChain(state_space, transition_matrix)

Creates a new discrete Markov chain object.

# Parameters
- `state_space`: The names of the states that make up the Markov chain.
- `transition_matrix`: The single step transition probability matrix.

# Examples
The following shows a basic Sunny-Cloudy-Rainy weather model.
```jldoctest core
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

```jldoctest core
println(transition_matrix(X))

# output

[0.9 0.1 0.0; 0.5 0.2 0.3; 0.1 0.4 0.5]
```

# References
https://en.wikipedia.org/wiki/Markov_chain#Discrete-time_Markov_chain
https://www.dartmouth.edu/~chance/teaching_aids/books_articles/probability_book/Chapter11.pdf
"""
struct DiscreteMarkovChain <: AbstractDiscreteMarkovChain
    state_space
    transition_matrix
    function DiscreteMarkovChain(state_space, transition_matrix)
        check(state_space, transition_matrix, DiscreteMarkovChain)
        new(state_space, transition_matrix)
    end
end
required_row_sum(::Core.Type{<:AbstractDiscreteMarkovChain}) = 1
function check(state_space, transition_matrix, type)
    if length(state_space) != size(transition_matrix)[1]
        error("The state space and transition matrix should be the same size.")
    end
    if length(unique(state_space)) != length(state_space)
        error("The state space must have unique elements.")
    end
    if !is_row_stochastic(transition_matrix, required_row_sum(type))
        error("The transition matrix should be row-stochastic 
        (each row must sum up to $(required_row_sum(type))).")
    end
end
function DiscreteMarkovChain(transition_matrix)
    return DiscreteMarkovChain(1:(size(transition_matrix)[1]), transition_matrix)
end

"""
    state_space(x)

# Parameters
- `x`: some kind of Markov chain.

# Returns
The state space of the Markov chain.
"""
state_space(x::AbstractMarkovChain) = x.state_space

"""
    transition_matrix(x)

# Parameters
- `x`: some kind of Markov chain.

# Returns
The transition matrix of the Markov chain.
"""
transition_matrix(x::AbstractMarkovChain) = x.transition_matrix

"""
    state_index(x)

# Parameters
- `x`: some kind of Markov chain.

# Returns 
A dictionary mapping each state in a Markov chain to its position in the state space. 
It is essentially the inverse of state_space(x).
"""
state_index(x::AbstractMarkovChain) = Dict(
    state => index for (index, state) in enumerate(state_space(x))
)
Base.length(x::AbstractMarkovChain) = length(state_space(x))

"""
    digraph(x)

Creates a digraph (directed graph) representation of a Markov chain.

# Parameters
- `x`: some kind of Markov chain.

# Returns
A 1D array of 2-tuples. An element ``(i, j)`` is in the array 
iff the transition matrix at ``i,j`` is nonzero.
"""
function digraph(x::AbstractMarkovChain)
    T = transition_matrix(x)
    V = 1:(size(T)[1])
    return [(i, j) for i in V for j in V if T[i, j] != 0] 
end

"""
    is_row_stochastic(mat, row_sum=1)

Tests whether a matrix, x, is row stochasitc.
The desired sum of each row can be specified as well.

# Parameters
- `mat`: a matrix that we want to check.
- `row_sum`: the desired value that each row should total to.

# Returns
`true` if the given matrix, `mat`, is row-stochasitc.
"""
function is_row_stochastic(mat, row_sum=1)
    n, p = size(mat)
    if n == 0
        return true
    end
    return repeat([row_sum], n) ≈ mat*ones(p)
end

"""
    communication_classes(x)

# Parameters
- `x`: some kind of Markov chain.

# Returns
A tuple containing 2 arrays. 
- This first array contains C arrays which store the states that communicate. 
- The second array is an array of Bool where the ith value is true if the 
  ith communication class is recurrent. 
"""
function communication_classes(x::AbstractMarkovChain)
    S = state_space(x)
    T = transition_matrix(x)

    V = 1:length(S)
    E = digraph(x)

    classes = strongly_connected_components(V, E)
    recurrence = []

    for class in classes
        submatrix = T[class, class]
        push!(recurrence, is_row_stochastic(submatrix, 1))
    end

    classes = [[S[i] for i in class] for class in classes]
    return classes, recurrence
end

"""
    periodicities(x::AbstractDiscreteMarkovChain)

A more advanced version of `communication_classes` 
designed for discrete Markov chains. It is the same as 
`communication_classes` but it returns periodicities as well.

# Parameters
- `x`: some kind of discrete Markov chain.

# Returns
A tuple containing 3 arrays. 
- This first array contains C arrays which store the states that communicate. 
- The second array is an array of Bool where the ith value is true if the 
  ith communication class is recurrent. 
- The third array is the periodicity of each communication class.
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
    decompose(x)

All transition matrices can be reordered so that we have
```math
T = \\begin{pmatrix}
A & 0\\\\
B & C
\\end{pmatrix}
```
Where ``A``, ``B`` and ``C`` are as described below. 
``0`` is a matrix of zeros.

# Parameters
- `x`: some kind of Markov chain.

# Returns 
A tuple of four values
- `states`: the first value is an array of the new states.
- `A`: the second value is a matrix of recurrent states to recurrent states.
- `B`: the third value is a matrix of transient states to recurrent states.
- `C`: the fourth value is a matrix of transient to transient states.
"""
function decompose(x::AbstractMarkovChain)
    T = transition_matrix(x)
    classes, recurrence = communication_classes(x)
    r_states = []
    t_states = []
    for (i, recurrent) in enumerate(recurrence)
        if recurrent
            append!(r_states, classes[i])
        else
            append!(t_states, classes[i])
        end
    end

    states = vcat(r_states, t_states)
    indexes = [state_index(x)[state] for state in states]

    A = [
        T[indexes[i], indexes[j]] 
        for i in 1:length(r_states), j in 1:length(r_states)
    ]
    B = [
        T[indexes[length(r_states)+i], indexes[j]] 
        for i in 1:length(t_states), j in 1:length(r_states)
    ]
    C = [
        T[indexes[length(r_states)+i], indexes[length(r_states)+j]] 
        for i in 1:length(t_states), j in 1:length(t_states)
    ]

    return states, A, B, C
end

"""
    canonical_form(x)

Reorders the states of the transition matrix of `x` so that 
recurrent states appear first and transient states appear last.


# Parameters
- `x`: some kind of Markov chain.

# Returns
A tuple with the new states and the new transition matrix.
"""
function canonical_form(x::AbstractMarkovChain)
    states, A, B, C = decompose(x)
    O = zeros(Int, size(A)[1], size(C)[2])
    result = vcat(
        hcat(A, O),
        hcat(B, C),
    )
    return states, result
end

"""
    is_regular(x)

# Parameters
- `x`: some kind of Markov chain.

# Returns
`true` if the Markov chain, `x`, is regular.
"""
function is_regular(x::AbstractMarkovChain)
    classes, _, periods = periodicities(x)
    if length(classes) == 0
        return false
    end
    return (length(classes) == 1) & (periods[1] == 1)
end

"""
    is_ergodic(x)

# Parameters
- `x`: some kind of Markov chain.

# Returns
`true` if the Markov chain, `x`, is ergodic.
This is, if every state can be accessed from every other state.
Another term for this is irreducible.
"""
function is_ergodic(x::AbstractMarkovChain)
    classes, _ = communication_classes(x)
    return length(classes) == 1
end

"""
    is_absorbing(x)

# Parameters
- `x`: some kind of Markov chain.

# Returns
`true` if the Markov chain, `x`, is an absorbing chain.
So the process is guarenteed to be absorbed eventually.
"""
function is_absorbing(x::AbstractMarkovChain)
    states, A, B, C = decompose(x)
    r = size(A)[1]
    return (r > 0) && (A ≈ LinearAlgebra.I)
end

"""
    stationary_distribution(x)

# Parameters
- `x`: some kind of Markov chain.

# Returns
A column vector, `w`, that satisfies the equation ``w'T = w``.
"""
function stationary_distribution(x::AbstractMarkovChain)
    T = transition_matrix(x)
    n = length(state_space(x))

    if n == 0
        return Any[]
    end

    a = (T - LinearAlgebra.I)'
    a[1, :] = ones(n)
    b = zeros(Int, n)
    b[1] = 1
    return a\b
end

"""
    fundamental_matrix(x)

The fundamental matrix of a markov chain is defined to be ``(I-C)^{-1}`` 
where ``C`` is the sub-transition matrix that takes transient states 
to transient states.

# Parameters
- `x`: some kind of Markov chain.

# Returns
The fundamental matrix of the Markov chain
"""
function fundamental_matrix(x::AbstractDiscreteMarkovChain)
    states, _, _, C = decompose(x)

    if size(C)[1] == 0
        return Array{Any}(undef, 0, 0)
    end

    return LinearAlgebra.inv(LinearAlgebra.I - C)
end

"""
    expected_time_to_absorption(x)

The expected number of steps that the process will take while in 
the transient super class given that the process started in state i.

# Parameters
- `x`: some kind of Markov chain.

# Returns
A 1D array where element ``i`` is the total number of revisits to 
transient state ``i`` before leaving the transient super class.

# Note
It is advised to have `x` in canonical form already.
This is to avoid confusion of what states make up each 
element of the array ouput.
"""
function expected_time_to_absorption(x::AbstractDiscreteMarkovChain)
    states, A, B, C = decompose(x)
    M = fundamental_matrix(x)
    EV = M*ones(Int, size(M)[1]) - ones(Int, size(M)[1])
    return EV
end

"""
    exit_probabilities(x)

# Parameters
- `x`: some kind of Markov chain.

# Returns
An array where element ``i, j`` is the probability that transient 
state ``i`` will enter recurrent state ``j`` on its first step 
out of the transient states. That is, ``e_{i,j}``.
"""
function exit_probabilities(x::AbstractDiscreteMarkovChain)
    M = fundamental_matrix(x)
    states, A, B, C = decompose(x)
    return M * B
end

"""
    first_passage_probabilities(x, t, i=missing, j=missing)

This is the probability that the process enters state ``j`` 
for the first time at time ``t`` given that the process started 
in state ``i`` at time 0. That is, ``f^{(t)}_{i,j}``. If no `i` 
or `j` is given, then it will return a matrix instead with 
entries ``f^{(t)}_{i,j}`` for `i` and `j` in the state space of `x`.

# Why Do We Use A Slow Algorithm?
So that `t` can be symbolic if nessesary. That is, if symbolic math
libraries want to use this library, it will pose no hassle.

# Parameters
- `x`: some kind of Markov chain.
- `t`: the time to calculate the first passage probability.
- `i`: the state that the prcess starts in.
- `j`: the state that the process must reach for the first time.

# Returns
A scalar value or a matrix depending on whether `i` and `j` are given.

# References
https://scholar.uwindsor.ca/cgi/viewcontent.cgi?article=1125&context=major-papers
http://maths.dur.ac.uk/stats/courses/ProbMC2H/_files/handouts/1516MarkovChains2H.pdf
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

    js = 1:n  # the columns to loop through
    calc_i_ne_j = true  # calculate the off-diagonals
    calc_i_eq_j = true  # calculate the diagonals
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

    Ft = zeros(eltype(T), n, n)  # empty matrix

    # if i != j
    if calc_i_ne_j
        for j in js

            P0 = copy(T)
            P0[1:n, j] = zeros(eltype(T), n, 1)
            F = P0^(t-1) * T
            Ft[1:n, j] = F[1:n, j]
        end
    end

    # if i == j
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

"""
    mean_recurrence_time(x)

This is the expected number of steps for the process to 
return to state i given that the process started in state i.

# Parameters
- `x`: some kind of Markov chain.

# Returns
A 1D array where the ith entry is the mean recurrence time of state i.

# Note
`x` must be irreducible (i.e. ergodic).
"""
function mean_recurrence_time(x::AbstractMarkovChain)
    return 1 ./ stationary_distribution(x)
end

"""
    mean_first_passage_time(x)

This is the expected number of steps for the process to 
reach state j given that the process started in state i.

# Parameters
- `x`: some kind of Markov chain.

# Returns
A matrix where the i,jth entry is the mean recurrence time of state i. 
Diagonal elements are 0. 
The diagonals would have represented the mean recurrence time.
"""
function mean_first_passage_time(x::AbstractMarkovChain)
    n = length(state_space(x))
    T = transition_matrix(x)

    if n == 0
        return T  # keeps eltype the same
    end

    W = repeat(stationary_distribution(x)', n)
    Z = LinearAlgebra.inv(LinearAlgebra.I - T + W)
    Zjj = repeat(LinearAlgebra.diag(Z)', n)

    M = (Zjj .- Z) ./ W
    return M
end

end # module
