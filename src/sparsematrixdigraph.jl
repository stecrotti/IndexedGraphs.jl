"""
    SparseMatrixDiGraph{T<:Integer} <: AbstractSparseMatrixGraph{T}
    
A type representing a sparse directed graph.

### FIELDS

- `A` -- square adjacency matrix filled with `NullNumber`s. Transpose of the adjacency matrix
- `X` -- square matrix for efficient access by row. `X[j,i]` points to the index of element `A[i,j]` in `A.nzval`. 
"""
struct SparseMatrixDiGraph{T<:Integer} <: AbstractSparseMatrixGraph{T}
    A :: SparseMatrixCSC{NullNumber, T}
    X :: SparseMatrixCSC{T, T}
end

function SparseMatrixDiGraph(A::SparseMatrixCSC)
    _checksquare(A)
    _check_selfloops(A)
    X = SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, collect(1:nnz(A)))
    At = sparse( permutedims(SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, fill(NullNumber(), length(A.nzval)))) )
    SparseMatrixDiGraph(At, X)
end

SparseMatrixDiGraph(A::AbstractMatrix) = SparseMatrixDiGraph(sparse(A))

Graphs.edges(g::SparseMatrixDiGraph) = (IndexedEdge{Int}(i, g.A.rowval[k], k) for i=1:size(g.A,2) for k=nzrange(g.A,i))

Graphs.ne(g::SparseMatrixDiGraph) = nnz(g.A)

# Graphs.outneighbors(g::SparseMatrixDiGraph, i::Integer) = @view g.A.rowval[nzrange(g.A,i)]
Graphs.inneighbors(g::SparseMatrixDiGraph, i::Integer) = @view g.X.rowval[nzrange(g.X,i)]

Graphs.is_directed(g::SparseMatrixDiGraph) = true

Graphs.is_directed(::Type{SparseMatrixDiGraph{T}}) where T = true

Base.zero(g::SparseMatrixDiGraph) = SparseMatrixDiGraph(zero(g.A))

outedges(g::SparseMatrixDiGraph, i::Integer) = (IndexedEdge{Int}(i, g.A.rowval[k], k) for k in nzrange(g.A, i))
inedges(g::SparseMatrixDiGraph, i::Integer) = (IndexedEdge{Int}(g.X.rowval[k], i, g.X.nzval[k]) for k in nzrange(g.X, i))

edge_idx(g::SparseMatrixDiGraph, src::Integer, dst::Integer) = nzindex(g.A, dst, src)
edge_src_dst(g::SparseMatrixDiGraph, id::Integer) = reverse(nzindex(g.A, id))

"""
    get_edge(g::SparseMatrixDiGraph, src::Integer, dst::Integer)
    get_edge(g::SparseMatrixDiGraph, id::Integer)

Get edge given source and destination or given edge index.
""" 
function get_edge(g::SparseMatrixDiGraph, src::Integer, dst::Integer)
    id = edge_idx(g, src, dst)
    IndexedEdge(src, dst, id)
end
function get_edge(g::SparseMatrixDiGraph, id::Integer)
    i, j = edge_src_dst(g, id)
    IndexedEdge(i, j, id)
end
