"""
    SparseMatrixDiGraph{T<:Integer} <: AbstractSparseMatrixGraph{T}
    
A type representing a sparse directed graph.

### FIELDS

- `A` -- square adjacency matrix filled with `NullNumber`s
- `X` -- square matrix for efficient access by row. `X[j,i]` points to the index of element `A[i,j]` in `A.nzval`. 
"""
struct SparseMatrixDiGraph{T<:Integer} <: AbstractSparseMatrixGraph{T}
    A :: SparseMatrixCSC{NullNumber, T}
    X :: SparseMatrixCSC{T, T}
end

function SparseMatrixDiGraph(A::AbstractMatrix{NullNumber})
    _checksquare(A)
    _check_selfloops(A)
    A = sparse(A)
    X = sparse(SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, collect(1:nnz(A)))')
    SparseMatrixDiGraph(A, X)
end

function SparseMatrixDiGraph(A::AbstractMatrix)
    A = sparse(A)
	SparseMatrixDiGraph(SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, fill(NullNumber(), length(A.nzval))))
end

Graphs.edges(g::SparseMatrixDiGraph) = (IndexedEdge{Int}(g.A.rowval[k], j, k) for j=1:size(g.A,2) for k=nzrange(g.A,j))

Graphs.ne(g::SparseMatrixDiGraph) = nnz(g.A)

Graphs.outneighbors(g::SparseMatrixDiGraph, i::Integer) = @view g.X.rowval[nzrange(g.X,i)]

Graphs.is_directed(g::SparseMatrixDiGraph) = true

Graphs.is_directed(::Type{SparseMatrixDiGraph{T}}) where T = true

Base.zero(g::SparseMatrixDiGraph) = SparseMatrixDiGraph(zero(g.A))

inedges(g::SparseMatrixDiGraph, i::Integer) = (IndexedEdge{Int}(g.A.rowval[k], i, k) for k in nzrange(g.A, i))

outedges(g::SparseMatrixDiGraph, i::Integer) = (IndexedEdge{Int}(i, g.X.rowval[k], g.X.nzval[k]) for k in nzrange(g.X, i))

edge_idx(g::SparseMatrixDiGraph, src::Integer, dst::Integer) = nzindex(g.A, src, dst)
edge_src_dst(g::SparseMatrixDiGraph, id::Integer) = nzindex(g.A, id)

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
