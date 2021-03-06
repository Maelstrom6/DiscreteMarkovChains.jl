module DiscreteMarkovChains
import LinearAlgebra

include("Utils.jl")

include("MC.jl")
include("DMC.jl")
include("CMC.jl")

export DiscreteMarkovChain, state_space, transition_matrix, probability_matrix,
communication_classes, periodicities, decompose, canonical_form,
is_regular, is_ergodic, is_absorbing, is_reversible,
stationary_distribution, fundamental_matrix,
mean_time_to_absorption, exit_probabilities, first_passage_probabilities,
mean_recurrence_time, mean_first_passage_time

export ContinuousMarkovChain, embedded

end  # Module
