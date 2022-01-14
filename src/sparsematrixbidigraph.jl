"""
    SparseMatrixBiDiGraph{T<:Integer} <: AbstractSparseMatrixDiGraph{T}
    
A type representing a sparse directed graph with access to both outedges and inedges.

### FIELDS

- `A` -- square matrix filled with `NullNumber`s. `A[i,j]` corresponds to an edge `j=>i`
- `X` -- square matrix for efficient access by row. `X[j,i]` points to the index of element `A[i,j]` in `A.nzval`. 
"""
struct SparseMatrixBiDiGraph{T<:Integer} <: AbstractSparseMatrixDiGraph{T}
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

Graphs.inneighbors(g::SparseMatrixBiDiGraph, i::Integer) = @view g.X.rowval[nzrange(g.X,i)]

Base.zero(g::SparseMatrixBiDiGraph) = SparseMatrixBiDiGraph(zero(g.A))

inedges(g::SparseMatrixBiDiGraph, i::Integer) = (IndexedEdge{Int}(g.X.rowval[k], i, g.X.nzval[k]) for k in nzrange(g.X, i))

# Returns sparse adj matrix. Elements default to Int (to match Graphs)
function Graphs.LinAlg.adjacency_matrix(g::SparseMatrixBiDiGraph, T::DataType=Int)
    SparseMatrixCSC(g.X.m, g.X.n, g.X.colptr, g.X.rowval, ones(T, nnz(g.X)))
end