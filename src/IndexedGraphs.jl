module IndexedGraphs

import SparseArrays:
    sparse, SparseMatrixCSC, nnz, nzrange, rowvals, nonzeros, spzeros

import Base: 
    show, ==

import Graphs:
    AbstractGraph, AbstractEdge, src, dst, edgetype, has_vertex, has_edge, ne, nv, 
    edges, vertices, neighbors, inneighbors, outneighbors, is_directed, is_bipartite,
    DijkstraState, dijkstra_shortest_paths

import Graphs.LinAlg:
    adjacency_matrix

import LinearAlgebra: 
    issymmetric

import TrackingHeaps: 
    TrackingHeap, pop!, NoTrainingWheels, MinHeapOrder

export
    idx, ==,
    AbstractIndexedGraph, inedges, outedges,
    # undirected graphs
    IndexedGraph, get_edge,
    # directed graphs
    AbstractIndexedDiGraph, IndexedDiGraph, IndexedBiDiGraph,
    # factor graphs
    FactorGraph, VariableOrFactor, Variable, Factor, FactorGraphEdge, 
    nvariables, nfactors, variables, factors, bipartite_view

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

src(e::IndexedEdge) = e.src
dst(e::IndexedEdge) = e.dst
idx(e::AbstractIndexedEdge) = e.idx

function Base.:(==)(e1::T, e2::T) where {T<:AbstractIndexedEdge}
    fns = fieldnames(T)
    all( getproperty(e1, fn) == getproperty(e2, fn) for fn in fns )
end

include("utils.jl")

include("abstractindexedgraph.jl")
include("indexedgraph.jl")
include("indexeddigraph.jl")
include("indexedbidigraph.jl")
include("factorgraph.jl")

include("algorithms/dijkstra.jl")

end # end module