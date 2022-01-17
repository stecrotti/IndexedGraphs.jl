using Documenter, IndexedGraphs
makedocs(sitename="IndexedGraphs Documentation", format = Documenter.HTML(prettyurls = false))

deploydocs(
    repo = "github.com/stecrotti/IndexedGraphs.jl.git",
)