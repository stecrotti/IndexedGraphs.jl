# IndexedGraph

A type representing undirected graphs.

```@docs
IndexedGraph
```

```@docs
IndexedGraph(A::AbstractMatrix)
```

```@docs
inedges(g::IndexedGraph, i::Integer)
outedges(g::IndexedGraph, i::Integer)
get_edge(g::IndexedGraph, src::Integer, dst::Integer)
```

## Overrides from Graphs.jl

```@docs
edges(g::IndexedGraph, i::Integer)
```