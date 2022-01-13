"""
    AbstractSparseMatrixGraph{T} <: AbstractGraph{T}
An abstract type representing a sparse graph structure.
`AbstractSparseMatrixGraph`s must have the following elements:
  - `A::SparseMatrixCSC` adjacency matrix
"""
abstract type AbstractSparseMatrixGraph{T} <: AbstractGraph{T} end

function Base.show(io::IO, g::AbstractSparseMatrixGraph)
    s = is_directed(g) ? "directed" : "undirected"
    println(io, "{$(nv(g)), $(ne(g))} ", s, " sparse ", eltype(g), " graph")
    Base.print_array(io, adjacency_matrix(g))
end

Base.eltype(g::AbstractSparseMatrixGraph{T}) where T = T

Graphs.edgetype(g::AbstractSparseMatrixGraph{T}) where T = IndexedEdge{T}

Graphs.has_vertex(g::AbstractSparseMatrixGraph, i::Integer) = i ≤ size(g.A, 2)

Graphs.has_edge(g::AbstractSparseMatrixGraph, i::Integer, j::Integer) = g.A[i,j] != 0

Graphs.nv(g::AbstractSparseMatrixGraph) = size(g.A, 2)

Graphs.vertices(g::AbstractSparseMatrixGraph) = 1:size(g.A, 2)

Graphs.inneighbors(g::AbstractSparseMatrixGraph, i::Integer) = @view g.A.rowval[nzrange(g.A,i)]

# Returns sparse adj matrix. Elements default to Int (to match Graphs)
function Graphs.LinAlg.adjacency_matrix(g::AbstractSparseMatrixGraph, T::DataType=Int)
    SparseMatrixCSC(g.A.m, g.A.n, g.A.colptr, g.A.rowval, ones(T, nnz(g.A)))
end