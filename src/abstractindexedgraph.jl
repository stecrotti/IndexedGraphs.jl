"""
    AbstractIndexedGraph{T} <: AbstractGraph{T}
An abstract type representing an indexed graph.
`AbstractIndexedGraph`s must have the following elements:
  - `A::SparseMatrixCSC` adjacency matrix
"""
abstract type AbstractIndexedGraph{T} <: AbstractGraph{T} end

function Base.show(io::IO, g::AbstractIndexedGraph)
    s = is_directed(g) ? "directed" : "undirected"
    println(io, "{$(nv(g)), $(ne(g))} ", s, " sparse ", eltype(g), " graph")
    Base.print_array(io, adjacency_matrix(g))
end

Base.eltype(g::AbstractIndexedGraph{T}) where T = T

edgetype(g::AbstractIndexedGraph{T}) where T = IndexedEdge{T}

has_vertex(g::AbstractIndexedGraph, i::Integer) = i ≤ size(g.A, 2)

has_edge(g::AbstractIndexedGraph, i::Integer, j::Integer) = g.A[i,j] != 0

nv(g::AbstractIndexedGraph) = size(g.A, 2)

vertices(g::AbstractIndexedGraph) = 1:size(g.A, 2)

outneighbors(g::AbstractIndexedGraph, i::Integer) = @view g.A.rowval[nzrange(g.A,i)]

# Returns sparse adj matrix. Elements default to Int (to match Graphs)
function adjacency_matrix(g::AbstractIndexedGraph, T::DataType=Int)
    SparseMatrixCSC(g.A.m, g.A.n, g.A.colptr, g.A.rowval, ones(T, nnz(g.A)))
end