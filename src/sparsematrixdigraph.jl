"""
    AbstractSparseMatrixDiGraph{T} <: AbstractSparseMatrixGraph{T}

Abstract type for representing directed graphs.
"""
abstract type AbstractSparseMatrixDiGraph{T} <: AbstractSparseMatrixGraph{T}; end

Graphs.edges(g::AbstractSparseMatrixDiGraph) = (IndexedEdge{Int}(i, g.A.rowval[k], k) for i=1:size(g.A,2) for k=nzrange(g.A,i))
Graphs.ne(g::AbstractSparseMatrixDiGraph) = nnz(g.A)

Graphs.is_directed(g::AbstractSparseMatrixDiGraph) = true
Graphs.is_directed(::Type{AbstractSparseMatrixDiGraph{T}}) where T = true


outedges(g::AbstractSparseMatrixDiGraph, i::Integer) = (IndexedEdge{Int}(i, g.A.rowval[k], k) for k in nzrange(g.A, i))

edge_idx(g::AbstractSparseMatrixDiGraph, src::Integer, dst::Integer) = nzindex(g.A, dst, src)
edge_src_dst(g::AbstractSparseMatrixDiGraph, id::Integer) = reverse(nzindex(g.A, id))

"""
    get_edge(g::AbstractSparseMatrixDiGraph, src::Integer, dst::Integer)
    get_edge(g::AbstractSparseMatrixDiGraph, id::Integer)

Get edge given source and destination or given edge index.
""" 
function get_edge(g::AbstractSparseMatrixDiGraph, src::Integer, dst::Integer)
    id = edge_idx(g, src, dst)
    IndexedEdge(src, dst, id)
end
function get_edge(g::AbstractSparseMatrixDiGraph, id::Integer)
    i, j = edge_src_dst(g, id)
    IndexedEdge(i, j, id)
end

# Returns sparse adj matrix. Elements default to Int (to match Graphs)
function Graphs.LinAlg.adjacency_matrix(g::AbstractSparseMatrixDiGraph, T::DataType=Int)
    M = sparse(transpose(g.A))
    SparseMatrixCSC(M.m, M.n, M.colptr, M.rowval, ones(T, nnz(M)))
end

"""
    SparseMatrixDiGraph{T<:Integer} <: AbstractSparseMatrixDiGraph{T}
    
A type representing a sparse directed graph with access only to outedges.

### FIELDS

- `A` -- square matrix filled with `NullNumber`s. `A[i,j]` corresponds to an edge `j=>i`
"""
struct SparseMatrixDiGraph{T<:Integer} <: AbstractSparseMatrixDiGraph{T}
    A :: SparseMatrixCSC{NullNumber, T}
    function SparseMatrixDiGraph(A::SparseMatrixCSC{NullNumber, T}) where {T<:Integer}
        _checksquare(A)
        _check_selfloops(A)
        new{T}(A)
    end
end

function SparseMatrixDiGraph(A::AbstractMatrix)
    B = convert(SparseMatrixCSC, transpose(A))
    At = SparseMatrixCSC(B.m, B.n, B.colptr, B.rowval, fill(NullNumber(), nnz(B)))
    SparseMatrixDiGraph(At)
end

# WARNING: very slow! Not recommended, if you need inneighbors, check out `SparseMatrixBiDiGraph`
# here only to comply with the requirements for subtyping `AbstractGraph` from Graphs.jl
function Graphs.inneighbors(g::SparseMatrixDiGraph, i::Integer) 
    X = sparse(transpose(g.A))
    @view X.rowval[nzrange(X,i)]
end 

Base.zero(g::SparseMatrixDiGraph) = SparseMatrixDiGraph(zero(g.A))