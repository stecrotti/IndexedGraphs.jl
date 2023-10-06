# BipartiteIndexedGraph

A graph is [bipartite](https://en.wikipedia.org/wiki/Bipartite_graph) if the set of vertices
can be partitioned into two blocks such that there is no edge between vertices of different blocks.
Here we adopt the notation of referring to these as the "left" and "right" block.

A `BipartiteIndexedGraph` behaves just like a bipartite, undirected `IndexedGraph`, with two differences:
- Vertices can be indexed as usual via their integer index (which is called here a [`linearindex`](@ref)), or via a [`BipartiteGraphVertex`](@ref), i.e. by specifying `Left` or `Right` and the integer index of that vertex within its block. The typical use case is that where one has two vectors storing vertex properties of the two blocks, possibly with different `eltype`s, and each with indices starting at one.
- The adjacency matrix of a bipartite graph (possibly after permutation of the vertex indices) is made of two specular rectangular submatrices `A`. Only one of these is stored, leading to a slight improvement in efficiency.

`BipartiteIndexedGraph`s use the same edge type `IndexedEdge` as the other `AbstractIndexedGraph`s, which stores source and destination as linear indices. To retrieve the in-block indices of the two vertices incident on an edge, use [`vertex_left`](@ref), [`vertex_right`](@ref).

```@docs
BipartiteIndexedGraph
```

```@docs
BipartiteIndexedGraph(A::AbstractMatrix)
BipartiteIndexedGraph(g::AbstractGraph)
nv_left
nv_right
Left
Right
LeftorRight
BipartiteGraphVertex
vertex(i::Integer, ::Type{<:LeftorRight})
vertex(g::BipartiteIndexedGraph, i::Integer)
linearindex
vertices_left
vertices_right
vertex_left
vertex_right
inedges(g::BipartiteIndexedGraph, i::Integer)
inedges(g::BipartiteIndexedGraph, l::BipartiteGraphVertex{Left})
outedges(g::BipartiteIndexedGraph, i::Integer)
outedges(g::BipartiteIndexedGraph, l::BipartiteGraphVertex{Left})
```

## Overrides from Graphs.jl

```@docs
inneighbors(g::BipartiteIndexedGraph, i::Integer)
inneighbors(g::BipartiteIndexedGraph, l::BipartiteGraphVertex{Left})
outneighbors(g::BipartiteIndexedGraph, i::Integer)
outneighbors(g::BipartiteIndexedGraph, v::BipartiteGraphVertex)
adjacency_matrix(g::BipartiteIndexedGraph, T::DataType=Int)
```