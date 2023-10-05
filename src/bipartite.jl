struct BipartiteIndexedGraph{T<:Integer} <: AbstractIndexedGraph{T}
    A :: SparseMatrixCSC{NullNumber, T}
    X :: SparseMatrixCSC{T, T}
end

function BipartiteIndexedGraph(A::AbstractMatrix{NullNumber})
    A = sparse(A)
    X = sparse(SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, collect(1:nnz(A)))')
    BipartiteIndexedGraph(A, X)
end

function BipartiteIndexedGraph(A::AbstractMatrix)
    A = sparse(A)
	BipartiteIndexedGraph(SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, fill(NullNumber(), length(A.nzval))))
end

function BipartiteIndexedGraph(g::AbstractGraph)
    is_directed(g) && throw(ArgumentError("Only an undirected graph can be converted into a `BipartiteIndexedGraph`"))
    is_bipartite(g) || throw(ArgumentError("Graph must be bipartite"))
    bm = bipartite_map(g)
    p = vcat([i for i in eachindex(bm) if bm[i]==1], [i for i in eachindex(bm) if bm[i]==2])
    pl = [i for i in eachindex(bm) if bm[i]==1]
    pr = [i for i in eachindex(bm) if bm[i]==2]
    nl = sum(isequal(1), bm)
    A = adjacency_matrix(g)[p,p][nl+1:end, 1:nl]
    A = adjacency_matrix(g)[pl, pr]
    return BipartiteIndexedGraph(A)
end

function Base.show(io::IO, g::BipartiteIndexedGraph{T}) where T
    nl = nv_left(g)
    nr = nv_right(g)
    ned = ne(g)
    println(io, "BipartiteIndexedGraph{$T} with $nl + $nr vertices and $ned edges")
end

nv_left(g::BipartiteIndexedGraph) = size(g.A, 1)
nv_right(g::BipartiteIndexedGraph) = size(g.A, 2)
Graphs.ne(g::BipartiteIndexedGraph) = nnz(g.A)
Graphs.nv(g::BipartiteIndexedGraph) = nv_left(g) + nv_right(g)

Graphs.is_directed(g::BipartiteIndexedGraph) = false
Graphs.is_bipartite(g::BipartiteIndexedGraph) = true

struct Left end
struct Right end
LeftorRight = Union{Left, Right}

struct BipartiteGraphVertex{LR<:LeftorRight,T<:Integer}
    i  :: T
    function BipartiteGraphVertex{LR}(i::T) where {LR<:LeftorRight,T<:Integer}
        i > 0 || throw(ArgumentError("Vertex index must be positive, got $i")) 
        new{LR,T}(i)
    end
end

vertex(i::Integer, ::Type{LR}) where LR<:LeftorRight = BipartiteGraphVertex{LR}(i)
linearindex(g::BipartiteIndexedGraph, r::BipartiteGraphVertex{Right}) = r.i + nv_left(g)
linearindex(::BipartiteIndexedGraph, l::BipartiteGraphVertex{Left}) = l.i
function linearindex(g::BipartiteIndexedGraph, i::Integer, ::Type{LR}) where LR<:LeftorRight
    return linearindex(g, vertex(i, LR))
end

function vertex(g::BipartiteIndexedGraph, i::Integer)
    has_vertex(g, i) || throw(ArgumentError("Index $i not in range of vertices."))
    nl = nv_left(g)
    if 1 ≤ i ≤ nl
        return vertex(i, Left)
    else
        return vertex(i - nl, Right)
    end
end

Graphs.edgetype(::BipartiteIndexedGraph{T}) where T = IndexedEdge{T}

function Graphs.edges(g::BipartiteIndexedGraph{T}) where T
    (IndexedEdge(linearindex(g, g.A.rowval[k], Left), linearindex(g, j, Right), k) 
        for j=1:size(g.A,2) for k=nzrange(g.A,j))
end

Graphs.vertices(g::BipartiteIndexedGraph) = 1:nv(g)
vertices_left(g::BipartiteGraphVertex) = 1:nv_left(g)
vertices_right(g::BipartiteGraphVertex) = nv_left(g)+1:nv(g)
Graphs.has_vertex(g::BipartiteIndexedGraph, i::Integer) = i ∈ vertices(g)

Graphs.has_edge(g::BipartiteIndexedGraph, e::IndexedEdge) = e ∈ edges(g)
_checkrightindex(g::BipartiteIndexedGraph, i::Integer) = nv_left(g) + 1 ≤ i ≤ nv(g)
function Graphs.has_edge(g::BipartiteIndexedGraph, s::Integer, d::Integer)
    (has_vertex(g, s) && has_vertex(g, d)) || return false
    l, rlin = extrema((s, d))
    _checkrightindex(g, rlin) || return false
    r = vertex(g, rlin).i
    return !iszero(g.A[l, r])
end


function Graphs.inneighbors(g::BipartiteIndexedGraph, l::BipartiteGraphVertex{Left})
    return (linearindex(g, r, Right) for r in @view g.X.rowval[nzrange(g.X, l.i)])
end
function Graphs.inneighbors(g::BipartiteIndexedGraph, r::BipartiteGraphVertex{Right})
    return @view g.A.rowval[nzrange(g.A, r.i)]
end
Graphs.outneighbors(g::BipartiteIndexedGraph, v::BipartiteGraphVertex) = inneighbors(g, v)

Graphs.inneighbors(g::BipartiteIndexedGraph, i::Integer) = inneighbors(g, vertex(g, i))
Graphs.outneighbors(g::BipartiteIndexedGraph, i::Integer) = inneighbors(g, i)

function inedges(g::BipartiteIndexedGraph, l::BipartiteGraphVertex{Left})
    return (
        IndexedEdge(linearindex(g, g.X.rowval[k], Right), linearindex(g, l), g.X.nzval[k])
            for k in nzrange(g.X, l.i)
    )
end
function inedges(g::BipartiteIndexedGraph, r::BipartiteGraphVertex{Right})
    return (
        IndexedEdge(linearindex(g, g.A.rowval[k], Left), linearindex(g, r), k)
            for k in nzrange(g.A, r.i)
    )
end
inedges(g::BipartiteIndexedGraph, i::Integer) = inedges(g, vertex(g, i))

function outedges(g::BipartiteIndexedGraph, l::BipartiteGraphVertex{Left})
    return (
        IndexedEdge(linearindex(g, l), linearindex(g, g.X.rowval[k], Right), g.X.nzval[k])
            for k in nzrange(g.X, l.i)
    )
end
function outedges(g::BipartiteIndexedGraph, r::BipartiteGraphVertex{Right})
    return (
        IndexedEdge(linearindex(g, r), linearindex(g, g.A.rowval[k], Left), k)
            for k in nzrange(g.A, r.i)
    )
end
outedges(g::BipartiteIndexedGraph, i::Integer) = outedges(g, vertex(g, i))