using Documenter, IndexedGraphs
makedocs(sitename="IndexedGraphs Documentation", 
    # format = Documenter.HTML(prettyurls = false),
    pages = [
        "Home" => "index.md",
        "Graph types" => [
            "graph.md", 
            "digraph.md",
            "bidigraph.md"            
        ]
    ]
    )

deploydocs(
    repo = "github.com/stecrotti/IndexedGraphs.jl.git",
)