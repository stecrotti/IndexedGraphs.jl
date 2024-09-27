"""
    BipartiteIndexedGraph{T<:Integer} <: AbstractIndexedGraph{T}
    
A type representing a sparse, undirected bipartite graph.

### FIELDS

- `A` -- adjacency matrix filled with `NullNumber`s. Rows are vertices belonging to the left block, columns to the right block
- `X` -- square matrix for efficient access by row. `X[j,i]` points to the index of element `A[i,j]` in `A.nzval`. 
"""
struct BipartiteIndexedGraph{T<:Integer} <: AbstractIndexedGraph{T}
    A :: SparseMatrixCSC{NullNumber, T}
    X :: SparseMatrixCSC{T, T}
end

function BipartiteIndexedGraph(A::AbstractMatrix{NullNumber})
    A = sparse(A)
    X = sparse(SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, collect(1:nnz(A)))')
    BipartiteIndexedGraph(A, X)
end

"""
    BipartiteIndexedGraph(A::AbstractMatrix)

Construct a `BipartiteIndexedGraph` from adjacency matrix `A` with the convention that
rows are vertices belonging to the left block, columns to the right block
"""
function BipartiteIndexedGraph(A::AbstractMatrix)
    A = sparse(A)
	BipartiteIndexedGraph(SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, fill(NullNumber(), length(A.nzval))))
end

"""
    BipartiteIndexedGraph(g::AbstractGraph)

Build a `BipartiteIndexedGraph` from any undirected, bipartite graph.
"""
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

"""
    nv_left(g::BipartiteIndexedGraph)

Return the number of vertices in the left block
"""
nv_left(g::BipartiteIndexedGraph) = size(g.A, 1)

"""
    nv_right(g::BipartiteIndexedGraph)

Return the number of vertices in the right block
"""
nv_right(g::BipartiteIndexedGraph) = size(g.A, 2)

Graphs.ne(g::BipartiteIndexedGraph) = nnz(g.A)
Graphs.nv(g::BipartiteIndexedGraph) = nv_left(g) + nv_right(g)

Graphs.is_directed(g::BipartiteIndexedGraph) = false
Graphs.is_directed(::Type{BipartiteIndexedGraph{T}}) where T = false
Graphs.is_bipartite(g::BipartiteIndexedGraph) = true

"""
    Left

Singleton type used to represent a vertex belonging to the left block in a [`BipartiteGraphVertex`](@ref)
"""
struct Left end

"""
    Right

Singleton type used to represent a vertex belonging to the ` Right` block in a [`BipartiteGraphVertex`](@ref)
"""
struct Right end

"""
    LeftorRight

`LeftorRight = Union{Left, Right}`
"""
LeftorRight = Union{Left, Right}

"""
    BipartiteGraphVertex

A `BipartiteGraphVertex{LR<:LeftorRight,T<:Integer}` represents a vertex in a bipartite graph.

### PARAMETERS
- `LR` -- Either [`Left`](@ref) or [`Right`](@ref)

### FIELDS
- `i` -- The index of the vertex within its block.
"""
struct BipartiteGraphVertex{LR<:LeftorRight,T<:Integer}
    i  :: T
    function BipartiteGraphVertex{LR}(i::T) where {LR<:LeftorRight,T<:Integer}
        i > 0 || throw(ArgumentError("Vertex index must be positive, got $i")) 
        new{LR,T}(i)
    end
end
"""
    vertex(i::Integer, ::Type{<:LeftorRight})

Build a [`BipartiteGraphVertex`](@ref)
"""
vertex(i::Integer, ::Type{LR}) where LR<:LeftorRight = BipartiteGraphVertex{LR}(i)

"""
    vertex(g::BipartiteIndexedGraph, i::Integer)

Build the [`BipartiteGraphVertex`](@ref) corresponding to linear index `i`.
Throws an error if `i` is not in the range of vertices of `g`
"""
function vertex(g::BipartiteIndexedGraph, i::Integer)
    has_vertex(g, i) || throw(ArgumentError("Index $i not in range of vertices."))
    nl = nv_left(g)
    if 1 ≤ i ≤ nl
        return vertex(i, Left)
    else
        return vertex(i - nl, Right)
    end
end

"""
    linearindex(g::BipartiteIndexedGraph, v::BipartiteGraphVertex{LR}) where {LR<:LeftorRight}
    linearindex(g::BipartiteIndexedGraph, i::Integer, ::Type{LR}) where LR<:LeftorRight

Return the linear index of a vertex, specified either by a [`BipartiteGraphVertex`](@ref) or by its index and block.
"""
linearindex(g::BipartiteIndexedGraph, r::BipartiteGraphVertex{Right}) = r.i + nv_left(g)
linearindex(::BipartiteIndexedGraph, l::BipartiteGraphVertex{Left}) = l.i
function linearindex(g::BipartiteIndexedGraph, i::Integer, ::Type{LR}) where LR<:LeftorRight
    return linearindex(g, vertex(i, LR))
end

Graphs.edgetype(::BipartiteIndexedGraph{T}) where T = IndexedEdge{T}

function Graphs.edges(g::BipartiteIndexedGraph{T}) where T
    (IndexedEdge(linearindex(g, g.A.rowval[k], Left), linearindex(g, j, Right), k) 
        for j=1:size(g.A,2) for k=nzrange(g.A,j))
end

Graphs.vertices(g::BipartiteIndexedGraph) = 1:nv(g)
Graphs.has_vertex(g::BipartiteIndexedGraph, i::Integer) = i ∈ vertices(g)

"""
    vertices_left(g::BipartiteIndexedGraph)

Return a lazy iterator to the vertices in the left block
"""
vertices_left(g::BipartiteIndexedGraph) = 1:nv_left(g)

"""
    vertices_right(g::BipartiteIndexedGraph)

Return a lazy iterator to the vertices in the right block
"""
vertices_right(g::BipartiteIndexedGraph) = nv_left(g)+1:nv(g)

_checkrightindex(g::BipartiteIndexedGraph, i::Integer) = nv_left(g) + 1 ≤ i ≤ nv(g)

Graphs.has_edge(g::BipartiteIndexedGraph, e::IndexedEdge) = e ∈ edges(g)
function Graphs.has_edge(g::BipartiteIndexedGraph, s::Integer, d::Integer)
    (has_vertex(g, s) && has_vertex(g, d)) || return false
    l, rlin = extrema((s, d))
    _checkrightindex(g, rlin) || return false
    r = vertex(g, rlin).i
    return !iszero(g.A[l, r])
end

_ordered_edge(e::IndexedEdge) = extrema((src(e), dst(e)))

"""
    vertex_left(g::BipartiteIndexedGraph, e::IndexedEdge)

Return the (in-block) index of the left-block vertex in `e`.
"""
function vertex_left(g::BipartiteIndexedGraph, e::IndexedEdge)
    l, rlin = _ordered_edge(e)
    return vertex(g, l).i
end

"""
    vertex_right(g::BipartiteIndexedGraph, e::IndexedEdge)

Return the (in-block) index of the right-block vertex in `e`.
"""
function vertex_right(g::BipartiteIndexedGraph, e::IndexedEdge)
    l, rlin = _ordered_edge(e)
    return vertex(g, rlin).i
end

""" 
    inneighbors(g::BipartiteIndexedGraph, v::BipartiteGraphVertex{<:LeftorRight})

Return a lazy iterator to the neighbors of variable `v` specified by a [`BipartiteGraphVertex`](@ref). 
"""
function Graphs.inneighbors(g::BipartiteIndexedGraph, l::BipartiteGraphVertex{Left})
    return (linearindex(g, r, Right) for r in @view g.X.rowval[nzrange(g.X, l.i)])
end
function Graphs.inneighbors(g::BipartiteIndexedGraph, r::BipartiteGraphVertex{Right})
    return @view g.A.rowval[nzrange(g.A, r.i)]
end

""" 
    inneighbors(g::BipartiteIndexedGraph, i::Integer)

Return a lazy iterator to the neighbors of variable `i` specified by its linear index. 
"""
Graphs.inneighbors(g::BipartiteIndexedGraph, i::Integer) = inneighbors(g, vertex(g, i))

""" 
    outneighbors(g::BipartiteIndexedGraph, v::BipartiteGraphVertex{<:LeftorRight})

Return a lazy iterator to the neighbors of variable `v` specified by a [`BipartiteGraphVertex`](@ref). 
"""
Graphs.outneighbors(g::BipartiteIndexedGraph, v::BipartiteGraphVertex) = inneighbors(g, v)

""" 
    outneighbors(g::BipartiteIndexedGraph, i::Integer)

Return a lazy iterator to the neighbors of variable `i` specified by its linear index. 
"""
Graphs.outneighbors(g::BipartiteIndexedGraph, i::Integer) = inneighbors(g, i)

"""
    degree(g::BipartiteIndexedGraph, v::BipartiteGraphVertex)

Return the degree of variable `v` specified by a [`BipartiteGraphVertex`](@ref). 
"""
Graphs.degree(g::BipartiteIndexedGraph, v::BipartiteGraphVertex) = length(inneighbors(g, v))

""" 
    inedges(g::BipartiteIndexedGraph, v::BipartiteGraphVertex{<:LeftorRight})

Return a lazy iterator to the edges incident on a variable `v` specified by a [`BipartiteGraphVertex`](@ref), with `v` as the destination. 
"""
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

""" 
    inedges(g::BipartiteIndexedGraph, i::Integer)

Return a lazy iterator to the edges incident on a variable `i` specified by its linear index, with `i` as the destination.
"""
inedges(g::BipartiteIndexedGraph, i::Integer) = inedges(g, vertex(g, i))


""" 
    outedges(g::BipartiteIndexedGraph, v::BipartiteGraphVertex{<:LeftorRight})

Return a lazy iterator to the edges incident on a variable `v` specified by a [`BipartiteGraphVertex`](@ref), with `v` as the source. 
"""
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

""" 
    outedges(g::BipartiteIndexedGraph, i::Integer)

Return a lazy iterator to the edges incident on a variable `i` specified by its linear index, with `i` as the source.
"""
outedges(g::BipartiteIndexedGraph, i::Integer) = outedges(g, vertex(g, i))

"""
    adjacency_matrix(g::BipartiteIndexedGraph, T::DataType=Int)

Return the symmetric adjacency matrix of size `nv(g) = nv_left(g) + nv_right(g)` 
where no distinction is made between left and right nodes.
"""
function Graphs.adjacency_matrix(g::BipartiteIndexedGraph, T::DataType=Int)
    m, n = nv_left(g), nv_right(g)
    A = SparseMatrixCSC(g.A.m, g.A.n, g.A.colptr, g.A.rowval, ones(T, nnz(g.A)))
    return [ spzeros(T, m, m)  A               ;
             A'                spzeros(T, n, n) ]
end


##### GENERATORS

"""
    rand_bipartite_graph([rng=default_rng()], nleft, nright, ned)

Create a bipartite graph with `nleft` nodes in the left block, `nright` nodes in the right block and `ned` edges taken uniformly at random.
"""
function rand_bipartite_graph(rng::AbstractRNG, nleft::Integer, nright::Integer, ned::Integer)
    nleft > 0 || throw(ArgumentError("Number of variable nodes must be positive, got $nleft"))
    nright > 0 || throw(ArgumentError("Number of left nodes must be positive, got $nright"))
    ned > 0 || throw(ArgumentError("Number of edges must be positive, got $ned"))
    nedmax = nleft * nright
    ned ≤ nedmax || throw(ArgumentError("Maximum number of edges is $nleft*$nright=$nedmax, got $ned"))

    I = zeros(Int, ned)
    J = zeros(Int, ned)
    K = ones(Int, ned)
    n = 1
    while n ≤ ned
        I[n] = rand(rng, 1:nleft)
        J[n] = rand(rng, 1:nright)
        if !any((i,j) == (I[n], J[n]) for (i,j) in Iterators.take(zip(I,J), n-1))
            n += 1
        end
    end
    A = sparse(I, J, K, nleft, nright)
    return BipartiteIndexedGraph(A)
end
function rand_bipartite_graph(nleft::Integer, nright::Integer, ned::Integer)
    rand_bipartite_graph(default_rng(), nleft, nright, ned)
end

"""
    rand_bipartite_graph([rng=default_rng()], nleft, nright, p)

Create a bipartite graph with `nleft` nodes in the left block, `nright` nodes in the right block and edges taken independently with probability `p`.
"""
function rand_bipartite_graph(rng::AbstractRNG, nright::Integer, nleft::Integer, p::Real)
    nright > 0 || throw(ArgumentError("Number of right nodes must be positive, got $nright"))
    nleft > 0 || throw(ArgumentError("Number of left nodes must be positive, got $nleft"))
    0 ≤ p ≤ 1 || throw(ArgumentError("Probability must be in [0,1], got $ned"))

    I = zeros(Int, 0)
    J = zeros(Int, 0)
    for (a, i) in Iterators.product(1:nleft, 1:nright)
        if rand(rng) < p
            push!(I, a)
            push!(J, i)
        end
    end
    K = ones(Int, length(I))
    A = sparse(I, J, K, nleft, nright)
    return BipartiteIndexedGraph(A)
end
function rand_bipartite_graph(nright::Integer, nleft::Integer, p::Real)
    rand_bipartite_graph(default_rng(), nleft, nright, p)
end

"""
    rand_regular_bipartite_graph([rng=default_rng()], nleft, nright, k)

Create a bipartite graph with `nleft` nodes in the left block, `nright` nodes in the right block, where all left nodes have degree `k`.
"""
function rand_regular_bipartite_graph(rng::AbstractRNG, nright::Integer, nleft::Integer, 
        k::Integer)
    nright > 0 || throw(ArgumentError("Number of right nodes must be positive, got $nright"))
    nleft > 0 || throw(ArgumentError("Number of left nodes must be positive, got $nleft"))
    k > 0 || throw(ArgumentError("Degree `k` must be positive, got $k"))
    k ≤ nright || throw(ArgumentError("Degree `k` must be smaller or equal than number of rights, got $k")) 

    I = reduce(vcat, fill(a, k) for a in 1:nleft)
    J = reduce(vcat, sample(rng, 1:nright, k; replace=false) for _ in 1:nleft)
    K = ones(Int, length(I))
    A = sparse(I, J, K, nleft, nright)
    return BipartiteIndexedGraph(A)
end
function rand_regular_bipartite_graph(nright::Integer, nleft::Integer, k::Integer)
    rand_regular_bipartite_graph(default_rng(), nleft, nright, k)
end

"""
    rand_bipartite_tree([rng=default_rng()], n)

Create a tree bipartite graph with `n` vertices in total. The proportion of left/right nodes is random.
"""
function rand_bipartite_tree(rng::AbstractRNG, n::Integer)
    gg = prufer_decode(rand(rng, 1:n, n-2))
    return BipartiteIndexedGraph(gg)
end
rand_bipartite_tree(n::Integer) = rand_bipartite_tree(default_rng(), n)