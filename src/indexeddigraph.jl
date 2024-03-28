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
function Graphs.inneighbors(g::IndexedDiGraph, i::Integer) 
    X = sparse(transpose(g.A))
    @view X.rowval[nzrange(X,i)]
end 

Base.zero(g::IndexedDiGraph) = IndexedDiGraph(zero(g.A))

function Base.show(io::IO, g::IndexedDiGraph{T}) where T
    println(io, "{$(nv(g)), $(ne(g))} IndexedDiGraph{$T}")
end
