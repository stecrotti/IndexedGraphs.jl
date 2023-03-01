"""
    AbstractIndexedDiGraph{T} <: AbstractIndexedGraph{T}

Abstract type for representing directed graphs.
"""
abstract type AbstractIndexedDiGraph{T} <: AbstractIndexedGraph{T}; end

edges(g::AbstractIndexedDiGraph) = @inbounds (IndexedEdge{Int}(i, g.A.rowval[k], k) for i=1:size(g.A,2) for k=nzrange(g.A,i))
ne(g::AbstractIndexedDiGraph) = nnz(g.A)

is_directed(g::AbstractIndexedDiGraph) = true
is_directed(::Type{AbstractIndexedDiGraph{T}}) where T = true

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
function adjacency_matrix(g::AbstractIndexedDiGraph, T::DataType=Int)
    M = sparse(transpose(g.A))
    SparseMatrixCSC(M.m, M.n, M.colptr, M.rowval, ones(T, nnz(M)))
end

"""
    IndexedDiGraph{T<:Integer} <: AbstractIndexedDiGraph{T}
    
A type representing a sparse directed graph with access only to outedges.

### FIELDS

- `A` -- square matrix filled with `NullNumber`s. `A[i,j]` corresponds to an edge `j=>i`
"""
struct IndexedDiGraph{T<:Integer} <: AbstractIndexedDiGraph{T}
    A :: SparseMatrixCSC{NullNumber, T}
    function IndexedDiGraph(A::SparseMatrixCSC{NullNumber, T}) where {T<:Integer}
        _checksquare(A)
        _check_selfloops(A)
        new{T}(A)
    end
end

"""
    IndexedDiGraph(A::AbstractMatrix)

Constructs a IndexedDiGraph from the adjacency matrix A.

`IndexedDiGraph` internally stores the transpose of `A`. To avoid overhead due
to the transposition, use `IndexedDiGraph(transpose(At))` where `At` is the 
transpose of `A`.
"""
function IndexedDiGraph(A::AbstractMatrix)
    B = convert(SparseMatrixCSC, transpose(A))
    At = SparseMatrixCSC(B.m, B.n, B.colptr, B.rowval, fill(NullNumber(), nnz(B)))
    IndexedDiGraph(At)
end

"""
    IndexedDiGraph(A::AbstractSimpleGraph)

Construct an `IndexedDiGraph` from any `AbstractSimpleGraph` (Graphs.jl), 
directed or otherwise.
"""
IndexedDiGraph(sg::AbstractSimpleGraph) = IndexedDiGraph(adjacency_matrix(sg))

# WARNING: very slow! Not recommended, if you need inneighbors, check out `IndexedBiDiGraph`
# here only to comply with the requirements for subtyping `AbstractGraph` from Graphs.jl
function inneighbors(g::IndexedDiGraph, i::Integer) 
    X = sparse(transpose(g.A))
    @view X.rowval[nzrange(X,i)]
end 

Base.zero(g::IndexedDiGraph) = IndexedDiGraph(zero(g.A))

function Base.show(io::IO, g::IndexedDiGraph{T}) where T
    println(io, "{$(nv(g)), $(ne(g))} IndexedDiGraph{$T}")
end