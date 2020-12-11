abstract type AbstractContinuousMarkovChain <: AbstractMarkovChain end

"""
    ContinuousMarkovChain(transition_matrix)
    ContinuousMarkovChain(state_space, transition_matrix)
    ContinuousMarkovChain(discrete_markov_chain)

Creates a new continuous Markov chain object.
This is also known as a Markov jump process.

Note that an irreducible finite continuous-time Markov chain is
always positive recurrent and its stationary distribution always
exists, is unique and is equal to the limiting distribution.

# Arguments
- `state_space`: The names of the states that make up the Markov chain.
- `transition_matrix`: The transition intensity matrix. Also known as the generator matrix.
- `discrete_markov_chain`: An instance of `DiscreteMarkovChain`.

# Examples
The following shows a basic Sunny-Cloudy-Rainy weather model.
```jldoctest ContinuousMarkovChain
using DiscreteMarkovChains
T = [
    -.1 0.1 0.0;
    0.5 -.8 0.3;
    0.1 0.4 -.5;
]
X = ContinuousMarkovChain(["Sunny", "Cloudy", "Rainy"], T)
println(state_space(X))

# output

["Sunny", "Cloudy", "Rainy"]
```

```jldoctest ContinuousMarkovChain
println(transition_matrix(X))

# output

[-0.1 0.1 0.0; 0.5 -0.8 0.3; 0.1 0.4 -0.5]
```

# References
1. [Wikipedia](https://en.wikipedia.org/wiki/Continuous-time_Markov_chain)
"""
struct ContinuousMarkovChain <: AbstractContinuousMarkovChain
    state_space
    transition_matrix
    function ContinuousMarkovChain(state_space, transition_matrix)
        check(state_space, transition_matrix, ContinuousMarkovChain)
        new(state_space, transition_matrix)
    end
end

required_row_sum(::Core.Type{<:AbstractContinuousMarkovChain}) = 0
function ContinuousMarkovChain(transition_matrix)
    return ContinuousMarkovChain(1:(size(transition_matrix)[1]), transition_matrix)
end

function characteristic_matrix(::AbstractContinuousMarkovChain)
    return LinearAlgebra.UniformScaling(0)
end

probability_matrix(x::AbstractContinuousMarkovChain) = exp(transition_matrix(x))

convert(
    ::Type{ContinuousMarkovChain},
    x::AbstractDiscreteMarkovChain,
) = ContinuousMarkovChain(x)
function ContinuousMarkovChain(x::AbstractDiscreteMarkovChain)
    S = state_space(x)
    T = transition_matrix(x)

    if length(S) == 0
        return ContinuousMarkovChain(S, T)
    end

    Q = log(T)
    return ContinuousMarkovChain(S, Q)
end

convert(
    ::Type{DiscreteMarkovChain},
    x::AbstractContinuousMarkovChain,
) = DiscreteMarkovChain(x)
function DiscreteMarkovChain(x::AbstractContinuousMarkovChain)
    S = state_space(x)
    Q = transition_matrix(x)

    if length(S) == 0
        return DiscreteMarkovChain(S, Q)
    end

    T = exp(Q)
    return DiscreteMarkovChain(S, T)
end

function embedded(x::AbstractContinuousMarkovChain)
    Q = transition_matrix(x)

    diag_Q_inv = LinearAlgebra.diagm(1 ./ LinearAlgebra.diag(Q))
    T = -diag_Q_inv*Q
    T[LinearAlgebra.diagind(T)] .= 0  # Sometimes NaNs appear

    return DiscreteMarkovChain(state_space(x), T)
end
