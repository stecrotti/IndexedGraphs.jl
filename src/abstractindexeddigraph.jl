"""
    AbstractIndexedDiGraph{T} <: AbstractIndexedGraph{T}

Abstract type for representing directed graphs.
"""
abstract type AbstractIndexedDiGraph{T} <: AbstractIndexedGraph{T}; end

Graphs.edges(g::AbstractIndexedDiGraph) = @inbounds (IndexedEdge{Int}(i, g.A.rowval[k], k) for i=1:size(g.A,2) for k=nzrange(g.A,i))
Graphs.ne(g::AbstractIndexedDiGraph) = nnz(g.A)

Graphs.is_directed(g::AbstractIndexedDiGraph) = true
Graphs.is_directed(::Type{<:AbstractIndexedDiGraph}) = true

"""
    outedges(g::AbstractIndexedDiGraph, i::Integer)

Return a lazy iterator to the edges outgoing from node `i` in `g`.
"""
outedges(g::AbstractIndexedDiGraph, i::Integer) = @inbounds (IndexedEdge{Int}(i, g.A.rowval[k], k) for k in nzrange(g.A, i))

edge_idx(g::AbstractIndexedDiGraph, src::Integer, dst::Integer) = nzindex(g.A, dst, src)
edge_src_dst(g::AbstractIndexedDiGraph, id::Integer) = reverse(nzindex(g.A, id))

"""
    get_edge(g::AbstractIndexedDiGraph, src::Integer, dst::Integer)
    get_edge(g::AbstractIndexedDiGraph, id::Integer)

Get edge given source and destination or given edge index.
""" 
function get_edge(g::AbstractIndexedDiGraph, src::Integer, dst::Integer)
    id = edge_idx(g, src, dst)
    IndexedEdge(src, dst, id)
end
function get_edge(g::AbstractIndexedDiGraph, id::Integer)
    i, j = edge_src_dst(g, id)
    IndexedEdge(i, j, id)
end

# Returns sparse adj matrix. Elements default to Int (to match Graphs)
function Graphs.LinAlg.adjacency_matrix(g::AbstractIndexedDiGraph, T::DataType=Int)
    M = sparse(transpose(g.A))
    SparseMatrixCSC(M.m, M.n, M.colptr, M.rowval, ones(T, nnz(M)))
end