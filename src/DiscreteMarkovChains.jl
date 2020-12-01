module DiscreteMarkovChains
import LinearAlgebra

include("Utils.jl")
export DiscreteMarkovChain, state_space, transition_matrix, 
communication_classes, periodicities, decompose, canonical_form,
is_regular, is_ergodic, is_absorbing,
stationary_distribution, fundamental_matrix, expected_time_to_absorption

abstract type AbstractMarkovChain end
abstract type AbstractDiscreteMarkovChain <: AbstractMarkovChain end

"""
    DiscreteMarkovChain(state_space, transition_matrix)

Creates a new discrete Markov chain object.

# Parameters
- `state_space`: The names of the states that make up the Markov chain.
- `transition_matrix`: The single step transition probability matrix.

# Examples
The following shows a basic Sunny-Cloudy-Rainy weather model.
```jldoctest
using DiscreteMarkovChains
T = [
     0.9 0.1 0;
     0.5 0.2 0.3;
     0.1 0.4 0.5
    ]
X = DiscreteMarkovChain(["Sunny", "Cloudy", "Rainy"], T)
print(state_space(X))

# output

["Sunny", "Cloudy", "Rainy"]
```
"""
struct DiscreteMarkovChain <: AbstractDiscreteMarkovChain
    state_space
    transition_matrix
end

"""
    state_space(x)

# Parameters
- `x`: some kind of Markov chain.

# Returns
The state space of a Markov chain.
"""
state_space(x::AbstractMarkovChain) = x.state_space

"""
    transition_matrix(x)

# Returns
The transition matrix of a Markov chain.
"""
transition_matrix(x::AbstractMarkovChain) = x.transition_matrix

"""
    state_index(x)

# Parameters
- `x`: some kind of Markov chain.

# Returns 
A dictionary mapping each state in a Markov chain to its position in the state space. It is essentially the inverse of state_space(x).
"""
state_index(x::AbstractMarkovChain) = Dict(state => index for (index, state) in enumerate(state_space(x)))
Base.length(x::AbstractMarkovChain) = length(state_space(x))

function digraph(x::AbstractMarkovChain)
    T = transition_matrix(x)
    V = 1:(size(T)[1])
    return [(i, j) for i in V for j in V if T[i, j] != 0] 
end

"""
    is_row_stochastic(x, row_sum=1)

Tests whether a matrix, x, is row stochasitc.
The desired sum of each row can be specified as well.

# Parameters
- `x`: some kind of Markov chain.

# Returns
`true` if the given matrix x is row-stochasitc.
"""
function is_row_stochastic(x, row_sum=1)
    n, p = size(x)
    return repeat([row_sum], n) ≈ x*ones(p)
end

"""
    communication_classes(x)

# Parameters
- `x`: some kind of Markov chain.

# Returns
A tuple containing 2 arrays. 
- This first array contains C arrays which store the states that communicate. 
- The second array is an array of Bool where the ith value is true if the ith communication class is recurrent. 
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
- The second array is an array of Bool where the ith value is true if the ith communication class is recurrent. 
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

# Parameters
- `x`: some kind of Markov chain.

# Returns 
A tuple of four values
- The first value is an array of the new states.
- The second value is a matrix of recurrent states to recurrent states.
- The third value is a matrix of transient states to recurrent states.
- The fourth value is a matrix of transient to transient states.
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

    A = [T[indexes[i], indexes[j]] for i in 1:length(r_states), j in 1:length(r_states)]
    B = [T[indexes[length(r_states)+i], indexes[j]] for i in 1:length(t_states), j in 1:length(r_states)]
    C = [T[indexes[length(r_states)+i], indexes[length(r_states)+j]] for i in 1:length(t_states), j in 1:length(t_states)]

    return states, A, B, C
end

"""
    canonical_form(x)

Reorders the states of the transition matrix of `x` so that we have recurrent states first and transient states last.

# Returns
A tuple with the new states and the new transition matrix.
"""
function canonical_form(x::AbstractMarkovChain)
    states, A, B, C = decompose(x)
    O = zeros(Int, size(A)[1], size(C)[2])
    result = vcat(
        hcat(A, O),
        hcat(B, C)
    )
    return states, result
end

"""
    is_regular(x)

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

# Returns
A column vector, `w`, that satisfies the equation ``w'T = w``.
"""
function stationary_distribution(x::AbstractMarkovChain)
    T = transition_matrix(x)
    n = length(state_space(x))

    if n == 0
        return Array{Any}(undef, 0, 0)
    end

    a = (T - LinearAlgebra.I)'
    a[1, :] = ones(n)
    b = zeros(Int, n)
    b[1] = 1
    return a\b
end

function fundamental_matrix(x::AbstractDiscreteMarkovChain)
    states, _, _, C = decompose(x)
    return LinearAlgebra.inv(C - LinearAlgebra.I)
end

function expected_time_to_absorption(x::AbstractDiscreteMarkovChain)
    states, A, B, C = decompose(x)
    M = fundamental_matrix(x)
    EV = M*ones(Int, shape(M)[1], 1) - ones(Int, shape(M)[1], 1)
    return EV
end

function exit_probabilities(x::AbstractDiscreteMarkovChain)
    M = fundamental_matrix(x)
    states, A, B, C = decompose(x)
    return M * B
end

function first_passage_probabilities(x::AbstractDiscreteMarkovChain, t, i=missing, j=missing)
    n = length(state_space(x))
    S = state_space(x)
    T = transition_matrix(x)

    if n == 0
        return transition_matrix(x)
    end

    js = 1:n  # the columns to loop through
        calc_i_ne_j = true  # calculate the off-diagonals
        calc_i_eq_j = true  # calculate the diagonals
        if (i !== missing) && (j !== missing)
            js = [S[j]]
            if i == j
                calc_i_ne_j = false
            else
                calc_i_eq_j = false
            end
        end

        Ft = zeros(Int, n, n)  # empty matrix

        # if i != j
        if calc_i_ne_j
            for j in js

                P0 = copy(T)
                P0[0:n, j] = zeros(Int, n, 1)
                F = P0^(t-1) * P
                Ft[1:n, j] = F[1:n, j]
            end
        end

        # if i == j
        if calc_i_eq_j
            for j in js

                P_ = copy(T)
                P_[j, 1:n] = zeros(Int, 1, n)

                Pnew = zeros(Int, 2*n, 2*n)
                Pnew[1:n, 1:n] = P
                Pnew[(n+1):(2*n), (n+1):(2*n)] = P_
                Pnew[n+j, 1:n] = P[j, 1:n]

                P0 = copy(T)
                P0[1:(2*n), j] = zeros(Int, 2*n, 1)

                F = P0^(t - 1) * Pnew

                Ft[j, j] = F[n+j, j]
            end
        end

        if (i !== missing) && (j !== missing)
            return Ft[i, j]
        end
        return Ft
end

# https://cran.r-project.org/web/packages/markovchain/vignettes/an_introduction_to_markovchain_package.pdf

# exit_probabilities
# first_passage_probabilities

# meanFirstPassageTime
# meanRecurrenceTime
# meanAbsorptionTime

end # module
