# FactorGraph

A type representing [Factor Graphs](https://en.wikipedia.org/wiki/Factor_graph).

```@docs
FactorGraph
```

```@docs
FactorGraph(A::AbstractMatrix)
```

```@docs
get_edge(g::FactorGraph)
bipartite_view(g::FactorGraph, T::DataType=Int)
```

## Overrides from Graphs.jl

```@docs
Graphs.neighbors(g::FactorGraph, v::Variable)
Graphs.neighbors(g::FactorGraph, f::Factor)
Graphs.edges(g::FactorGraph, v::Variable)
Graphs.edges(g::FactorGraph, f::Factor)
Graphs.LinAlg.adjacency_matrix(g::FactorGraph, T::DataType=Int)
```