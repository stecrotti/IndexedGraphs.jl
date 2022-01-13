"""
    SparseMatrixGraph{T<:Integer} <: AbstractSparseMatrixGraph{T}

A type representing a sparse undirected graph.

### FIELDS

- `A` -- square adjacency matrix. `A[i,j] == A[j,i]` contains the unique index associated to unidrected edge `(i,j)`
"""
struct SparseMatrixGraph{T<:Integer} <: AbstractSparseMatrixGraph{T}
    A :: SparseMatrixCSC{T, T}   # Adjacency matrix. Values are unique edge id's
    
    function SparseMatrixGraph(A::SparseMatrixCSC{T}) where T
        _checksquare(A)
        _check_selfloops(A)
        M = SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, zeros(Int, nnz(A)))
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

SparseMatrixGraph(A::AbstractMatrix) = SparseMatrixGraph(sparse(A))

function Graphs.edges(g::SparseMatrixGraph) 
    (IndexedEdge{Int}(i, g.A.rowval[k], g.A.nzval[k])
        for i=1:size(g.A,2) for k=nzrange(g.A,i) if i < g.A.rowval[k])
end

Graphs.ne(g::SparseMatrixGraph) = Int( nnz(g.A) / 2 ) 

Graphs.neighbors(g::SparseMatrixGraph, i::Integer) = Graphs.inneighbors(g, i)

Graphs.outneighbors(g::SparseMatrixGraph, i::Integer) = Graphs.inneighbors(g, i)

Graphs.is_directed(g::SparseMatrixGraph) = false

Graphs.is_directed(::Type{SparseMatrixGraph{T}}) where T = false

Base.zero(g::SparseMatrixGraph) = SparseMatrixGraph(zero(g.A))

function Graphs.edges(g::SparseMatrixGraph, i::Integer)
    (IndexedEdge(extrema((i, g.A.rowval[k]))..., g.A.nzval[k]) 
        for k in nzrange(g.A, i))
end

function outedges(g::SparseMatrixGraph, i::Integer)
    (IndexedEdge(i, g.A.rowval[k], g.A.nzval[k]) for k in nzrange(g.A, i))
end

function inedges(g::SparseMatrixGraph, i::Integer)
    (IndexedEdge(g.A.rowval[k], i, g.A.nzval[k]) for k in nzrange(g.A, i))
end

function edge_idx(g::SparseMatrixGraph, src::Integer, dst::Integer)
    k = nzindex(g.A, src, dst)
    g.A.nzval[k]
end
function edge_src_dst(g::SparseMatrixGraph, id::Integer)
    k = findfirst(isequal(id), g.A.nzval)
    i, j = nzindex(g.A, k)
    return extrema((i,j))    # return sorted        
end

"""
    get_edge(g::SparseMatrixGraph, src::Integer, dst::Integer)
    get_edge(g::SparseMatrixGraph, id::Integer)

Get edge given source and destination or given edge index.
""" 
function get_edge(g::SparseMatrixGraph, src::Integer, dst::Integer)
    id = edge_idx(g, src, dst)
    IndexedEdge(src, dst, id)
end

function get_edge(g::SparseMatrixGraph, id::Integer)
    i, j = edge_src_dst(g, id)
    IndexedEdge(i, j, id)
end


