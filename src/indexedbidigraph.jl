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
"""
function IndexedBiDiGraph(A::AbstractMatrix) 
    _checksquare(A)
    _check_selfloops(A)
    B = convert(SparseMatrixCSC, transpose(A))
    X = permutedims( SparseMatrixCSC(B.m, B.n, B.colptr, B.rowval, collect(1:nnz(B))) )
    At = SparseMatrixCSC(B.m, B.n, B.colptr, B.rowval, fill(NullNumber(), nnz(B)))
    IndexedBiDiGraph(At, X)
end

"""
    IndexedBiDiGraph(A::AbstractSimpleGraph)

Construct an `IndexedBiDiGraph` from any `AbstractSimpleGraph` (Graphs.jl), 
directed or otherwise.
"""
IndexedBiDiGraph(sg::AbstractSimpleGraph) = IndexedBiDiGraph(adjacency_matrix(sg))

function Base.show(io::IO, g::IndexedBiDiGraph{T}) where T
    println(io, "{$(nv(g)), $(ne(g))} IndexedBiDiGraph{$T}")
end

Graphs.inneighbors(g::IndexedBiDiGraph, i::Integer) = @inbounds @view g.X.rowval[nzrange(g.X,i)]

Base.zero(g::IndexedBiDiGraph) = IndexedBiDiGraph(zero(g.A))

"""
    inedges(g::AbstractIndexedBiDiGraph, i::Integer)

Return a lazy iterator to the edges ingoing to node `i` in `g`.
"""
inedges(g::IndexedBiDiGraph, i::Integer) = @inbounds (IndexedEdge{Int}(g.X.rowval[k], i, g.X.nzval[k]) for k in nzrange(g.X, i))

# Return a copy of the adjacency matrix with elements of type `T`
function Graphs.LinAlg.adjacency_matrix(g::IndexedBiDiGraph, T::DataType=Int)
    SparseMatrixCSC(g.X.m, g.X.n, g.X.colptr, g.X.rowval, ones(T, nnz(g.X)))
end

"""
    issymmetric(g::IndexedBiDiGraph) -> Bool

Test whether a directed graph is symmetric, i.e. for each directed edge `i=>j` there also exists the edge `j=>i`
"""
function LinearAlgebra.issymmetric(g::IndexedBiDiGraph)
    for i in vertices(g)
        ein = inedges(g, i)
        eout = outedges(g, i)
        length(ein) == length(eout) || return false
        for (ei, eo) in zip(ein, eout)
            src(ei) == dst(eo) || return false
        end
    end
    return true
end

"""
    bidirected_with_mappings(g::IndexedGraph) -> (gdir, dir2undir, undir2dir)

Construct an `IndexedBiDiGraph` `gdir` from an `IndexedGraph` `g` by building two directed edges per every undirected edge in `g`.
In addition, return two vectors containing mappings from the undirected edges of `g` to the corresponding directed edges of `gdir`.

### OUTPUT
- `gdir` -- The directed graph
- `dir2undir` -- A vector of integers mapping the indices of the directed edges of `gdir` to the corresponding undirected edges of `g`
- `undir2dir` -- A vector of vectors with two integers each mapping the indices of the undirected edges of `g` to the two corresponding directed edges of `gdir`

"""
function bidirected_with_mappings(g::IndexedGraph{T}) where {T<:Integer}
    gdir = IndexedBiDiGraph(g.A)
    dir2undir = zeros(T, ne(gdir))
    undir2dir = [zeros(T, 0) for _ in edges(g)]

    for i in vertices(gdir)
        for (dir, undir) in zip(inedges(gdir, i), inedges(g,i))
            dir2undir[idx(dir)] = idx(undir)
            push!(undir2dir[idx(undir)], idx(dir))
        end
    end
    return gdir, dir2undir, undir2dir
end