# DiscreteMarkovChains

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Maelstrom6.github.io/DiscreteMarkovChains.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Maelstrom6.github.io/DiscreteMarkovChains.jl/dev)
[![Build Status](https://github.com/Maelstrom6/DiscreteMarkovChains.jl/workflows/CI/badge.svg)](https://github.com/Maelstrom6/DiscreteMarkovChains.jl/actions)
[![Coverage](https://codecov.io/gh/Maelstrom6/DiscreteMarkovChains.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Maelstrom6/DiscreteMarkovChains.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

DiscreteMarkovChains is a package that supports various functions relating to discrete Markov chains. In particular, it deals with discrete-time discrete-space time-homogenous finite Markov chains.

This library also deals with continuous Markov chains. Any function in the documentation that takes "some kind of Markov chain" as an argument can be a `DiscreteMarkovChain` or a `ContinuousMarkovChain`. Sadly there are very few examples for continuous Markov chains but they operate in the same way as discrete Markov chains.

## Installation

`DiscreteMarkovChains` should be up on Julia's package registry.

Simply type `] add DiscreteMarkovChains` into the Julia REPL.

## Documentation

See [the documentation](https://Maelstrom6.github.io/DiscreteMarkovChains.jl/dev) hosted on GitHub Pages.

## Usage

### Discrete Time

We want to find out if this chain is an absorbing chain.

```julia
using DiscreteMarkovChains

transition_matrix = [
    0.0 1.0 0.0;
    0.5 0.0 0.5;
    0.0 1.0 0.0;
]
chain = DiscreteMarkovChain(transition_matrix)
is_absorbing(chain)

# output

false
```

Let's try find the communication classes, see if they are recurrent and what their periodicity is.

```julia
periodicities(chain)

# output

([[1, 2, 3]], Any[true], Any[2])
```

This means that we have one communication class with 3 recurrent states. Their periodicity is 2.

Since we have a single communication class, we can calculate the mean recurrence times.

```julia
mean_recurrence_time(chain)

# output

3-element Array{Float64,1}:
 4.0
 2.0
 4.0
```

So the first and third states take an average of 4 time steps to return to itself. The second state takes an average of 2 steps to return to itself.

### Continuous Time

There is support for continuous Markov chains as well.

```julia
generator = [
    -3 1 2;
    0 -1 1;
    1 1 -2;
]
chain = ContinuousMarkovChain(generator)

communication_classes(chain)

# output

([[1, 2, 3]], Any[true])
```

So we have one communication class that is recurrent.

Calculate the stationary distribution of the chain.

```julia
stationary_distribution(chain)

# output

3-element Array{Float64,1}:
 0.125
 0.5
 0.375
```

Calculate the mean first passage time of the chain.

```julia
round.(mean_first_passage_time(chain), digits=2)

# output

3×3 Array{Float64,2}:
 0.0  1.0  0.67
 3.0  0.0  1.0
 2.0  1.0  0.0
```

## Authors

- Chris du Plessis

## License

MIT
