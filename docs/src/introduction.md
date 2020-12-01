# Introduction

## Installation

`DiscreteMarkovChains` should be up on Julia's package registry.

Simply type

```julia
Import Pkg
Pkg.add("DiscreteMarkovChains")
```

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