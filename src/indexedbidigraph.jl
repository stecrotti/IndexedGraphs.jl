"""
    IndexedBiDiGraph{T<:Integer} <: AbstractIndexedDiGraph{T}
    
A type representing a sparse directed graph with access to both outedges and inedges.

### FIELDS

- `A` -- square matrix filled with `NullNumber`s. `A[i,j]` corresponds to edge `j=>i`.
- `X` -- square matrix for efficient access by row. `X[j,i]` points to the index of element `A[i,j]` in `A.nzval`. 
"""
struct IndexedBiDiGraph{T<:Integer} <: AbstractIndexedDiGraph{T}
    A :: SparseMatrixCSC{NullNumber, T}
    X :: SparseMatrixCSC{T, T}
end

"""
    IndexedBiDiGraph(A::AbstractMatrix) 

Construct an `IndexedBiDiGraph` from the adjacency matrix `A`. 

`IndexedBiDiGraph` internally stores the transpose of `A`. To avoid overhead due
to the transposition, use `IndexedBiDiGraph(transpose(At))` where `At` is the 
transpose of `A`.

```@example
using SparseArrays
At = sprand(100, 100, 0.1)           # At[i,j] corresponds to edge j=>i
g = IndexedBiDiGraph(transpose(At))  
g.A.rowval === At.rowval
```
"""
function IndexedBiDiGraph(A::AbstractMatrix) 
    _checksquare(A)
    _check_selfloops(A)
    B = convert(SparseMatrixCSC, transpose(A))
    X = permutedims( SparseMatrixCSC(B.m, B.n, B.colptr, B.rowval, collect(1:nnz(B))) )
    At = SparseMatrixCSC(B.m, B.n, B.colptr, B.rowval, fill(NullNumber(), nnz(B)))
    IndexedBiDiGraph(At, X)
end

Graphs.inneighbors(g::IndexedBiDiGraph, i::Integer) = @view g.X.rowval[nzrange(g.X,i)]

Base.zero(g::IndexedBiDiGraph) = IndexedBiDiGraph(zero(g.A))

"""
    inedges(g::AbstractIndexedDiGraph, i::Integer)

Return a lazy iterator to the edges ingoing to node `i` in `g`.
"""
inedges(g::IndexedBiDiGraph, i::Integer) = (IndexedEdge{Int}(g.X.rowval[k], i, g.X.nzval[k]) for k in nzrange(g.X, i))

# Return a copy of the adjacency matrix with elements of type `T`
function Graphs.LinAlg.adjacency_matrix(g::IndexedBiDiGraph, T::DataType=Int)
    SparseMatrixCSC(g.X.m, g.X.n, g.X.colptr, g.X.rowval, ones(T, nnz(g.X)))
end
