# Core Types

These are types that are most often used. They are simple and are used by all other functions.

There is basic support for continuous Markov chains. All applicable functions also work for continuous chains.

## Contents

```@contents
Pages = ["basic.md"]
```

## Index

```@index
Pages = ["basic.md"]
```

## Documentation

The main type is `DiscreteMarkovChain`:

```@docs
DiscreteMarkovChains.DiscreteMarkovChain
DiscreteMarkovChains.ContinuousMarkovChain

DiscreteMarkovChains.state_space
DiscreteMarkovChains.transition_matrix
DiscreteMarkovChains.probability_matrix

DiscreteMarkovChains.embedded
```
