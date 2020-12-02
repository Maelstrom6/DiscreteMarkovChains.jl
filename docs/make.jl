using DiscreteMarkovChains
using Documenter

makedocs(;
    modules=[DiscreteMarkovChains],
    authors="Chris du Plessis",
    repo="https://github.com/Maelstrom6/DiscreteMarkovChains.jl/blob/{commit}{path}#L{line}",
    sitename="DiscreteMarkovChains.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Maelstrom6.github.io/DiscreteMarkovChains.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "introduction.md",
        "Functions" => [
            "functions/basic.md",
            "functions/simple.md",
            "functions/communication.md",
            "functions/probability.md",
            "functions/mean.md",
            "functions/internal.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/Maelstrom6/DiscreteMarkovChains.jl",
)
