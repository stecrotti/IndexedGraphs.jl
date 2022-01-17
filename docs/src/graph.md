# IndexedGraph

A type representing undirected graphs.

```@docs
IndexedGraph
```

```@docs
IndexedGraph(A::AbstractMatrix)
```

```@docs
inedges
outedges
get_edge(g::IndexedGraph, src::Integer, dst::Integer)
```

## Overrides from Graphs.jl

```@docs
Graphs.edges(g::IndexedGraph, i::Integer)
```