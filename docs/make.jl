using GasphaseReactions
using Documenter

DocMeta.setdocmeta!(GasphaseReactions, :DocTestSetup, :(using GasphaseReactions); recursive=true)

makedocs(;
    modules=[GasphaseReactions],
    authors="Vinod Janardhanan",
    repo="https://github.com/vinodjanardhanan/GasphaseReactions.jl/blob/{commit}{path}#{line}",
    sitename="GasphaseReactions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vinodjanardhanan.github.io/GasphaseReactions.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => "gchem.md",
        "Execution" => "run.md"        
    ],
)

deploydocs(;
    repo="github.com/vinodjanardhanan/GasphaseReactions.jl",
    devbranch="main",
)
