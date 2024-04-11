module IndexedGraphs

using SparseArrays: sparse, SparseMatrixCSC, nnz, nzrange, rowvals, nonzeros, spzeros

using Graphs: Graphs, AbstractGraph, SimpleGraph, AbstractSimpleGraph, AbstractEdge, 
    src, dst, edgetype, has_vertex, has_edge, ne, nv,
    edges, vertices, neighbors, inneighbors, outneighbors, is_directed, is_bipartite,
    bipartite_map,  DijkstraState, dijkstra_shortest_paths
    
using Graphs.LinAlg

using LinearAlgebra: LinearAlgebra, issymmetric

using TrackingHeaps: TrackingHeap, pop!, NoTrainingWheels, MinHeapOrder

export
    # Graphs.jl
    src, dst, edgetype, has_vertex, has_edge, ne, nv, adjacency_matrix, degree,
    edges, vertices, neighbors, inneighbors, outneighbors, is_directed, is_bipartite,
    # Base
    ==, iterate,
    # IndexedGraphs
    AbstractIndexedGraph, inedges, outedges, idx,
    # undirected graphs
    IndexedGraph, get_edge,
    # directed graphs
    AbstractIndexedDiGraph, IndexedDiGraph, IndexedBiDiGraph,
    # bipartite graphs
    BipartiteIndexedGraph, Left, Right, LeftorRight, BipartiteGraphVertex,
    nv_left, nv_right, vertex, linearindex, vertices_left, vertices_right,
    vertex_left, vertex_right, 
    is_directed, issymmetric,
    bidirected_with_mappings

"""
    AbstractIndexedEdge{T<:Integer} <: AbstractEdge{T}

Abstract type for indexed edge.
`AbstractIndexedEdge{T}`s must have the following elements:
- `idx::T` integer positive index 
"""
abstract type AbstractIndexedEdge{T<:Integer} <: AbstractEdge{T}; end
    
"""
    IndexedEdge{T<:Integer} <: AbstractIndexedEdge{T}

Edge type for `IndexedGraph`s. Edge indices can be used to access edge 
properties stored in separate containers.
"""
struct IndexedEdge{T<:Integer} <: AbstractIndexedEdge{T}
	src::T
	dst::T
	idx::T
end

Graphs.src(e::IndexedEdge) = e.src
Graphs.dst(e::IndexedEdge) = e.dst
idx(e::AbstractIndexedEdge) = e.idx

function Base.:(==)(e1::T, e2::T) where {T<:AbstractIndexedEdge}
    fns = fieldnames(T)
    all( getproperty(e1, fn) == getproperty(e2, fn) for fn in fns )
end

function Base.show(io::IO, e::AbstractIndexedEdge)
    print(io, "Indexed Edge $(src(e)) => $(dst(e)) with index $(idx(e))")
end

Base.iterate(e::IndexedEdge, args...) = iterate((e.src, e.dst, e.idx), args...)

include("utils.jl")
include("abstractindexedgraph.jl")
include("indexedgraph.jl")
include("abstractindexeddigraph.jl")
include("indexeddigraph.jl")
include("indexedbidigraph.jl")
include("bipartiteindexedgraph.jl")

include("algorithms/dijkstra.jl")

end # end module
