# DiscreteMarkovChains

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Maelstrom6.github.io/DiscreteMarkovChains.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Maelstrom6.github.io/DiscreteMarkovChains.jl/dev)
[![Build Status](https://github.com/Maelstrom6/DiscreteMarkovChains.jl/workflows/CI/badge.svg)](https://github.com/Maelstrom6/DiscreteMarkovChains.jl/actions)
[![Coverage](https://codecov.io/gh/Maelstrom6/DiscreteMarkovChains.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Maelstrom6/DiscreteMarkovChains.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

DiscreteMarkovChains is a package that supports various functions relating to discrete Markov chains. In particular, it deals with discrete-time discrete-space time-homogenous finite Markov chains.

A 100% pure-Julia stack library for functions and queries related to discrete Markov chains.

## Installation

`DiscreteMarkovChains` should be up on Julia's package registry.

Simply type `] add DiscreteMarkovChains` into the Julia REPL.

## Documentation

See [the documentation](https://Maelstrom6.github.io/DiscreteMarkovChains.jl/dev) hosted on GitHub Pages.

## Usage

We want to find out if this chain is an absorbing chain.

```jldoctest home; output = false
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

```jldoctest home; output = false
periodicities(chain)

# output

([[1, 2, 3]], Any[true], Any[2])
```

This means that we have one communication class with 3 recurrent states. Their periodicity is 2.

Since we have a single communication class, we can calculate the mean recurrence times.

```jldoctest home; output = false
mean_recurrence_time(chain)

# output

3-element Array{Float64,1}:
 4.0
 2.0
 4.0
```

So the first and third states take an average of 4 time steps to return to itself. The second state takes an average of 2 steps to return to itself.

## Authors

- Chris du Plessis

## License

MIT
