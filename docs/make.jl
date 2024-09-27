using Documenter, IndexedGraphs, Graphs, Random

makedocs(sitename="IndexedGraphs Documentation", 
    # format = Documenter.HTML(prettyurls = false),
    pages = [
        "Home" => "index.md",
        "Graph types" => [
            "graph.md", 
            "bidigraph.md",
            "digraph.md",
            "bipartite.md"          
        ],
        "Reference" => "reference.md"
    ]
    )

deploydocs(
    repo = "github.com/stecrotti/IndexedGraphs.jl.git",
    push_preview = true,
)