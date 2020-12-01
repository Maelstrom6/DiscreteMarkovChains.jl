# Functions

```@index
```

## Core Types

The main type is `DiscreteMarkovChain`:

```@docs
DiscreteMarkovChains.DiscreteMarkovChain
DiscreteMarkovChains.state_space
DiscreteMarkovChains.transition_matrix
```

## Simple Functions



```@docs
DiscreteMarkovChains.is_regular
DiscreteMarkovChains.is_ergodic
DiscreteMarkovChains.is_absorbing
```

## Communicaton Functions

```@docs
DiscreteMarkovChains.communication_classes
DiscreteMarkovChains.periodicities
DiscreteMarkovChains.decompose
DiscreteMarkovChains.canonical_form
```

## Probability Functions

```@docs
DiscreteMarkovChains.stationary_distribution
DiscreteMarkovChains.fundamental_matrix
DiscreteMarkovChains.expected_time_to_absorption
DiscreteMarkovChains.exit_probabilities
DiscreteMarkovChains.first_passage_probabilities
```
