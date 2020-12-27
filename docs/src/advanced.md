# Advanced

## Custom Element Types

This library also supports transition matrices whose elements are user-defined objects or even matrices. See `test/Advanced.jl` for examples.

## Conversion Between Discrete And Continuous

There is conversion ability between types. This can be done via `convert` or simply calling the class to convert to:

```jldoctest conversion
using DiscreteMarkovChains

T = [
    0.3 0.7 0.0;
    0.5 0.0 0.5;
    0.0 1.0 0.0;
]
X = DiscreteMarkovChain(T)

probability_matrix(X)

# output

3×3 Array{Float64,2}:
 0.3  0.7  0.0
 0.5  0.0  0.5
 0.0  1.0  0.0
```

```jldoctest conversion
Y = convert(ContinuousMarkovChain, X)
Y = ContinuousMarkovChain(X)  # Same operation

round.(abs.(probability_matrix(Y)), digits=2)

# output

3×3 Array{Float64,2}:
 0.3  0.7  0.0
 0.5  0.0  0.5
 0.0  1.0  0.0
```

The conversion is done so that the one-step transition probability matrix remains the same. If you want the embedded chain, then simply use `embedded(Y)`.
