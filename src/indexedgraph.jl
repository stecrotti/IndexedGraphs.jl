"""
    IndexedGraph{T<:Integer} <: AbstractIndexedGraph{T}

A type representing a sparse undirected graph.

### FIELDS

- `A` -- square adjacency matrix. `A[i,j] == A[j,i]` contains the unique index associated to unidrected edge `(i,j)`
"""
struct IndexedGraph{T<:Integer} <: AbstractIndexedGraph{T}
    A :: SparseMatrixCSC{T, T}   # Adjacency matrix. Values are unique edge id's

    function IndexedGraph(A::SparseMatrixCSC{T}) where T
        _checksquare(A)
        _check_selfloops(A)
        M = SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, ones(Int, nnz(A)))
        idx_map = inverse_edge_idx(M)
        cnt = 1
        # give increasing indices to upper triangular elements
        for j in 1:size(M, 2)
            for k in nzrange(M, j)
                i = M.rowval[k]
                if i < j
                    M.nzval[k] = cnt; cnt += 1
                else
                    break
                end
            end
        end
        @assert cnt - 1 == (nnz(M) / 2)  # sanity check
        # give (i=>j) the same idx as (j=>i)
        for j in 1:size(M, 2)
            for k in nzrange(M, j)
                i = M.rowval[k]
                if i > j
                    M.nzval[k] = M.nzval[idx_map[k]]
                end
            end
        end
        new{Int}(M)
    end
end

"""
    IndexedGraph(A::AbstractMatrix)

Construct an `IndexedGraph` from symmetric adjacency matrix A.
"""
IndexedGraph(A::AbstractMatrix) = IndexedGraph(convert(SparseMatrixCSC, A))

function edges(g::IndexedGraph) 
    (IndexedEdge{Int}(extrema((i, g.A.rowval[k]))..., g.A.nzval[k])
        for i=1:size(g.A,2) for k=nzrange(g.A,i) if i > g.A.rowval[k])
end

ne(g::IndexedGraph) = Int( nnz(g.A) / 2 ) 

neighbors(g::IndexedGraph, i::Integer) = outneighbors(g, i)

inneighbors(g::IndexedGraph, i::Integer) = outneighbors(g, i)

is_directed(g::IndexedGraph) = false

is_directed(::Type{IndexedGraph{T}}) where T = false

Base.zero(g::IndexedGraph) = IndexedGraph(zero(g.A))

"""
    edges(g::IndexedGraph, i::Integer)

Return a lazy iterators to the edges incident to `i`.

By default unordered edges sort source and destination nodes in increasing order.
See [outedges](@ref) and [inedges](@ref) if you need otherwise.
"""
function edges(g::IndexedGraph, i::Integer)
    (IndexedEdge(extrema((i, g.A.rowval[k]))..., g.A.nzval[k]) 
        for k in nzrange(g.A, i))
end

"""
    outedges(g::IndexedGraph, i::Integer)

Return a lazy iterators to the edges incident to `i` with `i` as the source.
"""
function outedges(g::IndexedGraph, i::Integer)
    (IndexedEdge(i, g.A.rowval[k], g.A.nzval[k]) for k in nzrange(g.A, i))
end

"""
    inedges(g::IndexedGraph, i::Integer)

Return a lazy iterators to the edges incident to `i` with `i` as the destination.
"""
function inedges(g::IndexedGraph, i::Integer)
    (IndexedEdge(g.A.rowval[k], i, g.A.nzval[k]) for k in nzrange(g.A, i))
end

function edge_idx(g::IndexedGraph, src::Integer, dst::Integer)
    k = nzindex(g.A, src, dst)
    g.A.nzval[k]
end
function edge_src_dst(g::IndexedGraph, id::Integer)
    k = findfirst(isequal(id), g.A.nzval)
    i, j = nzindex(g.A, k)
    return extrema((i,j))    # return sorted        
end

"""
    get_edge(g::IndexedGraph, src::Integer, dst::Integer)
    get_edge(g::IndexedGraph, id::Integer)

Get edge given source and destination or given edge index.
""" 
function get_edge(g::IndexedGraph, src::Integer, dst::Integer)
    id = edge_idx(g, src, dst)
    IndexedEdge(src, dst, id)
end

function get_edge(g::IndexedGraph, id::Integer)
    i, j = edge_src_dst(g, id)
    IndexedEdge(i, j, id)
end


