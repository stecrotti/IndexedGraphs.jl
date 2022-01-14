"""
    SparseMatrixBiDiGraph{T<:Integer} <: AbstractSparseMatrixGraph{T}
    
A type representing a sparse directed graph.

### FIELDS

- `A` -- square matrix filled with `NullNumber`s. `A[i,j]` corresponds to an edge `j=>i`
- `X` -- square matrix for efficient access by row. `X[j,i]` points to the index of element `A[i,j]` in `A.nzval`. 
"""
struct SparseMatrixBiDiGraph{T<:Integer} <: AbstractSparseMatrixGraph{T}
    A :: SparseMatrixCSC{NullNumber, T}
    X :: SparseMatrixCSC{T, T}
end

function SparseMatrixBiDiGraph(A::AbstractMatrix) 
    _checksquare(A)
    _check_selfloops(A)
    B = convert(SparseMatrixCSC, transpose(A))
    X = permutedims( SparseMatrixCSC(B.m, B.n, B.colptr, B.rowval, collect(1:nnz(B))) )
    At = SparseMatrixCSC(B.m, B.n, B.colptr, B.rowval, fill(NullNumber(), nnz(B)))
    SparseMatrixBiDiGraph(At, X)
end

Graphs.edges(g::SparseMatrixBiDiGraph) = (IndexedEdge{Int}(i, g.A.rowval[k], k) for i=1:size(g.A,2) for k=nzrange(g.A,i))

Graphs.ne(g::SparseMatrixBiDiGraph) = nnz(g.A)

# Graphs.outneighbors(g::SparseMatrixBiDiGraph, i::Integer) = @view g.A.rowval[nzrange(g.A,i)]
Graphs.inneighbors(g::SparseMatrixBiDiGraph, i::Integer) = @view g.X.rowval[nzrange(g.X,i)]

Graphs.is_directed(g::SparseMatrixBiDiGraph) = true

Graphs.is_directed(::Type{SparseMatrixBiDiGraph{T}}) where T = true

Base.zero(g::SparseMatrixBiDiGraph) = SparseMatrixBiDiGraph(zero(g.A))

outedges(g::SparseMatrixBiDiGraph, i::Integer) = (IndexedEdge{Int}(i, g.A.rowval[k], k) for k in nzrange(g.A, i))
inedges(g::SparseMatrixBiDiGraph, i::Integer) = (IndexedEdge{Int}(g.X.rowval[k], i, g.X.nzval[k]) for k in nzrange(g.X, i))

edge_idx(g::SparseMatrixBiDiGraph, src::Integer, dst::Integer) = nzindex(g.A, dst, src)
edge_src_dst(g::SparseMatrixBiDiGraph, id::Integer) = reverse(nzindex(g.A, id))

"""
    get_edge(g::SparseMatrixBiDiGraph, src::Integer, dst::Integer)
    get_edge(g::SparseMatrixBiDiGraph, id::Integer)

Get edge given source and destination or given edge index.
""" 
function get_edge(g::SparseMatrixBiDiGraph, src::Integer, dst::Integer)
    id = edge_idx(g, src, dst)
    IndexedEdge(src, dst, id)
end
function get_edge(g::SparseMatrixBiDiGraph, id::Integer)
    i, j = edge_src_dst(g, id)
    IndexedEdge(i, j, id)
end

# Returns sparse adj matrix. Elements default to Int (to match Graphs)
function Graphs.LinAlg.adjacency_matrix(g::SparseMatrixBiDiGraph, T::DataType=Int)
    SparseMatrixCSC(g.X.m, g.X.n, g.X.colptr, g.X.rowval, ones(T, nnz(g.X)))
end