# FactorGraph

A type representing [Factor Graphs](https://en.wikipedia.org/wiki/Factor_graph).

```@docs
FactorGraph
```

```@docs
FactorGraph(A::AbstractMatrix)
```

```@meta
using IndexedGraphs.VariableOrFactor
```

```@docs
get_edge(g::FactorGraph, x::VariableOrFactor, y::VariableOrFactor)
bipartite_view(g::FactorGraph, T::DataType=Int)
```

## Overrides from Graphs.jl

```@docs
neighbors(g::FactorGraph, v::Variable)
neighbors(g::FactorGraph, f::Factor)
edges(g::FactorGraph, v::Variable)
edges(g::FactorGraph, f::Factor)
adjacency_matrix(g::FactorGraph, T::DataType=Int)
```