using Documenter
using CartesianRuns

makedocs(
    sitename = "CartesianRuns.jl",
    modules  = [CartesianRuns],
    pages    = ["Home" => "index.md"],
    checkdocs = :exports,
)

deploydocs(
    repo      = "github.com/vlc1/CartesianRuns.jl",
    devbranch = "main",
)
