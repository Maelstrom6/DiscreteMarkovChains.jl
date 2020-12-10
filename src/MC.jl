abstract type AbstractMarkovChain end

function check(state_space, transition_matrix, type)
    if length(state_space) != size(transition_matrix)[1]
        error(
            "The state space, $(state_space), and "*
            "transition matrix should be the same size."
        )
    end
    if length(unique(state_space)) != length(state_space)
        error("The state space, $(state_space), must have unique elements.")
    end
    if !is_row_stochastic(transition_matrix, required_row_sum(type))
        error(
            "The transition matrix, $(transition_matrix), should be row-stochastic "*
            "(each row must sum up to $(required_row_sum(type)))."
        )
    end
end

"""
    state_space(x)

# Definitions
The state space of a Markov chain is the (ordered)
set of values that the process is able to take on.
For example, in a Sunny-Cloudy-Rainy weather model,
the state space is `["Sunny", "Cloudy", "Rainy"]`.

# Arguments
- `x`: some kind of Markov chain.

# Returns
The state space of the Markov chain.
"""
state_space(x::AbstractMarkovChain) = x.state_space

"""
    transition_matrix(x)

# Definitions
The one-step transition matrix, ``T``, of a discrete Markov chain
or the generator matrix of a continuous Markov chain.

# Arguments
- `x`: some kind of Markov chain.

# Returns
The transition matrix of the Markov chain.
"""
transition_matrix(x::AbstractMarkovChain) = x.transition_matrix

"""
    probability_matrix(x)

# Definitions
The one-step transition matrix, ``T``, of a Markov chain, ``\\{X_t\\}``
is a matrix whose ``(i,j)``th entry is the probability of the process
being in state ``j`` at time 1 given that the process started
in state ``i`` at time 0. That is
```math
T = \\{\\{p^{(1)}_{i,j}\\}\\}_{i,j ∈ S} = \\{\\{\\mathbb{P}(X_1=j | X_0=i)\\}\\}_{i,j ∈ S}
```
where ``S`` is the state space of ``\\{X_t\\}``.

# Arguments
- `x`: some kind of Markov chain.

# Returns
The one-step probability matrix of the Markov chain.
"""
probability_matrix(x::AbstractMarkovChain) = x.transition_matrix

"""
    state_index(x)

# Arguments
- `x`: some kind of Markov chain.

# Returns
A dictionary mapping each state in a Markov chain to its position in the state space.
It is essentially the inverse of `state_space(x)`.
"""
state_index(x::AbstractMarkovChain) = Dict(
    state => index for (index, state) in enumerate(state_space(x))
)

"""
    digraph(x)

Creates a digraph (directed graph) representation of a Markov chain.

# Arguments
- `x`: some kind of Markov chain.

# Returns
A 1D array of 2-tuples. An element ``(i, j)`` is in the array
iff the transition matrix at ``(i,j)`` is nonzero.
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

# Definitions
A matrix is said to be row stochasic if all its rows sum to 1.
This definition is extened so that all its rows sum to `row_sum`.

# Arguments
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
    # We only want to throw when we are sure it is incorrect.
    # Float16 seemed like a good compromise.
    isapprox(repeat([row_sum], n), mat*ones(p), atol=eps(Float16))
end

"""
    communication_classes(x)

# Definitions
A state ``j`` is accessible from state ``i`` if it is possible
to eventually reach state ``j`` given that the process started
in state ``i``. That is, ``∃ \\, t ∈ \\mathbb{N}`` such that
``p^{(t)}_{i,j} > 0``.

States ``i`` and ``j`` communicate if ``i`` is accessible
from ``j`` and ``j`` is accessible from ``i``.

A communication class is the set of states in a Markov chain
that are all accessible from one another. Communication classes
form a class in the mathematical sense. They also form a
partition of the state space.

A state ``i`` is recurrent if the process is guarenteed to
eventually return to state ``i`` given that the process
started in state ``i``. If a state ``i`` is not recurrent,
then it is said to be transient.

# Arguments
- `x`: some kind of Markov chain.

# Returns
A tuple containing 2 arrays.
- This first array contains C arrays which store the states that communicate.
- The second array is an array of Bool where the ith value is true if the
  ith communication class is recurrent.

# Examples
```jldoctest communication_classes
using DiscreteMarkovChains
T = [
    1 0;
    0 1;
]
X = DiscreteMarkovChain(T)

communication_classes(X)

# output

([[1], [2]], Any[true, true])
```

So the Markov chain has two communication classes and both are recurrent.

```jldoctest communication_classes
T = [
    .5 .5 .0;
    .2 .8 .0;
    .1 .2 .7;
]
X = DiscreteMarkovChain(["Sunny", "Cloudy", "Rainy"], T)

communication_classes(X)

# output

([["Sunny", "Cloudy"], ["Rainy"]], Any[true, false])
```

So the Sunny and Cloudy states communicate and are recurrent.
The Rainy state is transient.
Note that this is not a very good weather model since once it
stops raining, it will never rain again.
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
        push!(recurrence, is_row_stochastic(submatrix, required_row_sum(typeof(x))))
    end

    classes = [[S[i] for i in class] for class in classes]
    return classes, recurrence
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

# Arguments
- `x`: some kind of Markov chain.

# Returns
A tuple of four values
- `states`: the first value is an array of the new states.
- `A`: the second value is a matrix of recurrent states to recurrent states.
- `B`: the third value is a matrix of transient states to recurrent states.
- `C`: the fourth value is a matrix of transient to transient states.

# Examples
```jldoctest decompose
using DiscreteMarkovChains
T = [
    0.0 1.0 0.0;
    1.0 0.0 0.0;
    0.1 0.2 0.7;
]
X = DiscreteMarkovChain(["Sunny", "Cloudy", "Rainy"], T)

decompose(X)

# output

(Any["Sunny", "Cloudy", "Rainy"], [0.0 1.0; 1.0 0.0], [0.1 0.2], [0.7])
```
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

# Arguments
- `x`: some kind of Markov chain.

# Returns
A tuple with the new states and the new transition matrix.

# Examples
We will use the same matrix as previous examples but change the order of it.

```jldoctest canonical_form
using DiscreteMarkovChains
T = [
    0.7 0.2 0.1;
    0.0 0.0 1.0;
    0.0 1.0 0.0;
]
X = DiscreteMarkovChain(["Rainy", "Cloudy", "Sunny"], T)

canonical_form(X)

# output

(Any["Cloudy", "Sunny", "Rainy"], [0.0 1.0 0.0; 1.0 0.0 0.0; 0.2 0.1 0.7])
```

Here, Cloudy and Sunny are recurrent states so they appear first.
Rainy is a transient state so it appears last.
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
    is_ergodic(x)

# Definitions
A Markov chain is called an ergodic chain if it is possible to go
from every state to every state (not necessarily in one move).
That is, every state communicates and the whole tranistion
matrix forms one communication class.

In many books, ergodic Markov chains are called irreducible.

If a Markov chain is not irreducible then it is reducible.
This means that the Markov chain can be broken down into
smaller irreducible chains that do not communicate. One
chain might still be accessible from another though.

# Arguments
- `x`: some kind of Markov chain.

# Returns
`true` if the Markov chain, `x`, is ergodic.
This is, if every state can be accessed from every other state.
Another term for this is irreducible.

# Examples
We will set up a matrix with 2 communication classes and
show that it is not ergodic.

```jldoctest is_ergodic
using DiscreteMarkovChains
T = [
    0 1 0;
    1 0 0;
    0 0 1;
]
X = DiscreteMarkovChain(T)

is_ergodic(X)

# output

false
```

Repeat the above but now all states communicate.

```jldoctest is_ergodic
T = [
    0.0 0.5 0.5;
    0.0 0.0 1.0;
    1.0 0.0 0.0;
]
X = DiscreteMarkovChain(T)

is_ergodic(X)

# output

true
```

Notice how a periodic chain is regular no matter the periodicity.

```jldoctest is_ergodic
T = [
    0 1 0;
    0 0 1;
    1 0 0;
]
X = DiscreteMarkovChain(T)

is_ergodic(X)

# output

true
```
"""
function is_ergodic(x::AbstractMarkovChain)
    classes, _ = communication_classes(x)
    return length(classes) == 1
end

"""
    is_absorbing(x)

# Definitions
A Markov chain is absorbing if it has at least one absorbing
state, and if from every state it is possible to go to an absorbing
state (not necessarily in one step).

# Arguments
- `x`: some kind of Markov chain.

# Returns
`true` if the Markov chain, `x`, is an absorbing chain.
So the process is guarenteed to be absorbed eventually.

# Examples
The following is a typical example of an absorbing chain.

```jldoctest is_absorbing
using DiscreteMarkovChains
T = [
    1.0 0.0 0.0;
    0.1 0.8 0.1;
    0.0 0.3 0.7;
]
X = DiscreteMarkovChain(T)

is_absorbing(X)

# output

true
```

If a Markov chain does not have an absorbing state then the chain
is not absorbing. The converse is not true as we will see soon.

```jldoctest is_absorbing
T = [
    0.5 0.5 0.0;
    0.5 0.0 0.5;
    0.0 0.5 0.5;
]
X = DiscreteMarkovChain(T)

is_absorbing(X)

# output

false
```

In the following, the chain has multiple absorbing states but the
process is not guarenteed to be absorbed into those states.

```jldoctest is_absorbing
T = [
    1 0 0 0;
    0 1 0 0;
    0 0 0 1;
    0 0 1 0;
]
X = DiscreteMarkovChain(T)

is_absorbing(X)

# output

false
```

# References
1. [Dartmouth College](https://www.dartmouth.edu/~chance/teaching_aids/books_articles/probability_book/Chapter11.pdf)
"""
function is_absorbing(x::AbstractMarkovChain)
    states, A, B, C = decompose(x)
    r = size(A)[1]
    return (r > 0) && (LinearAlgebra.diag(A) ≈ repeat([required_row_sum(typeof(x))], r))
end

"""
    stationary_distribution(x)

# Definitions
A stationary distribution of a Markov chain is a probability distribution
that remains unchanged in the Markov chain as time progresses.
It is a row vector, ``w`` such that its elements sum to 1 and it satisfies
``wT = w``. ``T`` is the one-step transiton matrix of the Markov chain.

In other words, ``w`` is invariant by the matrix ``T``.

For simplicity, this function returns a column vector instead of a row vector.

For continuous Markov chains, the stationary_distribution is given by the
solution to ``wQ = 0`` where ``Q`` is the transition intensity matrix.

# Arguments
- `x`: some kind of Markov chain.

# Returns
A column vector, ``w``, that satisfies the equation ``w'T = w'``.

# Examples
The stationary distribution will always exist. However, it might not be unique.

If it is unique there are no problems.
```jldoctest stationary_distribution
using DiscreteMarkovChains
T = [
    4 2 4;
    1 0 9;
    3 5 2;
]//10
X = DiscreteMarkovChain(T)

stationary_distribution(X)

# output

3-element Array{Rational{Int64},1}:
 35//129
 12//43
 58//129
```

If there are infinite solutions then the principle solution is taken
(every free variable is set to 0). A Moore-Penrose inverse is used.

```jldoctest stationary_distribution
T = [
    0.4 0.6 0.0;
    0.6 0.4 0.0;
    0.0 0.0 1.0;
]
X = DiscreteMarkovChain(T)

stationary_distribution(X)

# output

3-element Array{Float64,1}:
 0.33333333333333337
 0.33333333333333337
 0.33333333333333337
```

# References
1. [Brilliant.org](https://brilliant.org/wiki/stationary-distributions/#:~:text=A%20stationary%20distribution%20of%20a,transition%20matrix%20P%2C%20it%20satisfies)
"""
function stationary_distribution(x::AbstractMarkovChain)
    T = transition_matrix(x)
    n = length(state_space(x))

    if n == 0
        return Any[]
    end

    a = (T - characteristic_matrix(x))'
    a[1, :] = ones(n)
    b = zeros(Int, n)
    b[1] = 1
    try
        return a\b
    catch e
        # This is a retarded way to handle this specific
        # exception but I don't know how else to deal with it.
        if e isa LinearAlgebra.SingularException
            # The pinv returns floats even if inputs are
            # rational. We will thus use it as a last resort.
            return LinearAlgebra.pinv(a)*b
        else
            rethrow(e)
        end
    end
end

"""
    mean_recurrence_time(x)

# Definitions
This is the expected number of steps for the process to
return to state i given that the process started in state ``i``.

# Arguments
- `x`: some kind of Markov chain.

# Returns
A 1D array where the ith entry is the mean recurrence time of state ``i``.

# Notes
`x` must be irreducible (i.e. ergodic).

# Examples
```jldoctest mean_recurrence_time
using DiscreteMarkovChains
T = [
    0.1 0.2 0.7;
    0.3 0.0 0.7;
    0.4 0.4 0.2;
]
X = DiscreteMarkovChain(T)

mean_recurrence_time(X)

# output

3-element Array{Float64,1}:
 3.4615384615384603
 4.090909090909091
 2.1428571428571432
```

If we have a reducible chain, then no error will be
raised but the output will be nonsense.

```jldoctest mean_recurrence_time
T = [
    1.0 0.0 0.0;
    0.1 0.5 0.4;
    0.0 0.3 0.7;
]
X = DiscreteMarkovChain(T)

mean_recurrence_time(X)

# output

3-element Array{Float64,1}:
   1.0
 -Inf
 -Inf
```
"""
function mean_recurrence_time(x::AbstractMarkovChain)
    return 1 ./ stationary_distribution(x)
end

"""
    mean_first_passage_time(x)

# Definitions
This is the expected number of steps for the process to
reach state ``j`` given that the process started in state ``i``.
We say that the mean first passage time between a state and itself is 0.

# Arguments
- `x`: some kind of Markov chain.

# Returns
A matrix where the ``(i,j)``th entry is the mean recurrence time of state ``i``.
Diagonal elements are 0.
The diagonals would have represented the mean recurrence time.

# Notes
`x` must be irreducible (i.e. ergodic).

# Examples
```jldoctest mean_first_passage_time
using DiscreteMarkovChains
T = [
    0.1 0.2 0.7;
    0.3 0.0 0.7;
    0.4 0.4 0.2;
]
X = DiscreteMarkovChain(T)

mean_first_passage_time(X)

# output

3×3 Array{Float64,2}:
 0.0      3.40909  1.42857
 2.88462  0.0      1.42857
 2.69231  2.95455  0.0
```

If we have a reducible chain, then no error will be
raised but the output will be nonsense.

```jldoctest mean_first_passage_time
T = [
    1.0 0.0 0.0;
    0.1 0.5 0.4;
    0.0 0.3 0.7;
]
X = DiscreteMarkovChain(T)

mean_first_passage_time(X)

# output

3×3 Array{Float64,2}:
  0.0     -Inf  -Inf
 23.3333  NaN   -Inf
 26.6667  -Inf  NaN
```
"""
function mean_first_passage_time(x::AbstractMarkovChain)
    n = length(state_space(x))
    T = probability_matrix(x)

    if n == 0
        return T  # Keeps eltype the same
    end

    W = repeat(stationary_distribution(x)', n)
    Z = LinearAlgebra.inv(LinearAlgebra.I - T + W)
    Zjj = repeat(LinearAlgebra.diag(Z)', n)

    M = (Zjj .- Z) ./ W
    return M
end
