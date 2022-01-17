"""
    FactorGraph{T<:Integer} <: AbstractIndexedGraph{T}
    
A type representing a sparse factor graph.

### FIELDS

- `A` -- adjacency matrix filled with `NullNumber`s. Rows are factors, columns are variables
- `X` -- square matrix for efficient access by row. `X[j,i]` points to the index of element `A[i,j]` in `A.nzval`. 
"""
struct FactorGraph{T<:Integer} <: AbstractIndexedGraph{T}
    A :: SparseMatrixCSC{NullNumber, T}
    X :: SparseMatrixCSC{T, T}
end

function FactorGraph(A::AbstractMatrix{NullNumber})
    _check_selfloops(A)
    A = sparse(A)
    X = sparse(SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, collect(1:nnz(A)))')
    FactorGraph(A, X)
end

"""
    FactorGraph(A::AbstractMatrix)

Construct a factor graph from adjacency matrix `A` with the convention that
rows are factors, columns are variables
"""
function FactorGraph(A::AbstractMatrix)
    A = sparse(A)
	FactorGraph(SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, fill(NullNumber(), length(A.nzval))))
end

abstract type VariableOrFactor{T<:Integer}; end

"""
    Variable{T<:Integer}

Wraps an index to specify that it is the index of a variable. 
See e.g. [`Graphs.neighbors(g::FactorGraph, v::Variable)`](@ref)
"""
struct Variable{T<:Integer} <: VariableOrFactor{T}
    i::T
    function Variable(i::T) where {T<:Integer}
        i > 0 || throw(ArgumentError("Variable index must be positive, got $i")) 
        new{T}(i)
    end
end

"""
    Factor{T<:Integer}

Wraps an index to specify that it is the index of a factor. 
See e.g. [`Graphs.neighbors(g::FactorGraph, f::Factor)`](@ref)
"""
struct Factor{T<:Integer} <: VariableOrFactor{T}
    a::T
    function Factor(a::T) where {T<:Integer}
        a > 0 || throw(ArgumentError("Factor index must be positive, got $a")) 
        new{T}(a)
    end
end

"""
FactorGraphEdge{T<:Integer} <: AbstractIndexedEdge{T}

Edge type for `FactorGraph`s.
### Safe constructors
    FactorGraphEdge(f::Factor, v::Variable, idx::Integer)
    FactorGraphEdge(v::Variable, f::Factor, idx::Integer)

Example:
```@example
    FactorGraphEdge(Factor(2), Variable(4), 3)
```
"""
struct FactorGraphEdge{T<:Integer} <: AbstractIndexedEdge{T}
	a::T    # factor index
	i::T    # variable index
    idx::T  # edge index
    function FactorGraphEdge(a::T, i::T, idx::T) where {T<:Integer}
        i > 0 || throw(ArgumentError("Variable index must be positive, got $i")) 
        a > 0 || throw(ArgumentError("Factor index must be positive, got $a")) 
        idx > 0 || throw(ArgumentError("Edge index must be positive, got $idx")) 
        new{T}(a, i, idx)
    end
end

function FactorGraphEdge(a::Integer, i::Integer, idx::Integer)
    a_, i_, idx_ = promote(a, i, idx)
    FactorGraphEdge(a_, i_, idx_)
end

function FactorGraphEdge(f::Factor{T}, v::Variable{T}, idx::T) where {T<:Integer}
    FactorGraphEdge(f.a, v.i, idx)
end
function FactorGraphEdge(v::Variable{T}, f::Factor{T}, idx::T) where {T<:Integer}
    FactorGraphEdge(f, v, idx)
end

function Base.show(io::IO, e::FactorGraphEdge)
    print(io, "Edge from factor $(e.a) to variable $(e.i) with index $(e.idx)")
end

Base.eltype(g::FactorGraph{T}) where T = T
Graphs.edgetype(g::FactorGraph{T}) where T = FactorGraphEdge{T}

nvariables(g::FactorGraph) = size(g.A, 2)
nfactors(g::FactorGraph) = size(g.A, 1)

factors(g::FactorGraph) = (Factor(a) for a in 1:nfactors(g))
variables(g::FactorGraph) = (Variable(i) for i in 1:nvariables(g))

Graphs.ne(g::FactorGraph) = nnz(g.A)
Graphs.nv(g::FactorGraph) = nvariables(g) + nfactors(g)

Graphs.edges(g::FactorGraph) = (FactorGraphEdge(g.A.rowval[k], j, k) for j=1:size(g.A,2) for k=nzrange(g.A,j))
# Graphs.vertices(g::FactorGraph) = Iterators.flatten((factors(g), variables(g)))
Graphs.vertices(g::FactorGraph) = 1:nv(g)

"""
    Graphs.edges(g::FactorGraph, v::Variable)

Return a lazy iterator to the edges incident to variable `v`
"""
function Graphs.edges(g::FactorGraph, v::Variable)
    i = v.i
    (FactorGraphEdge( g.A.rowval[k], i, k) for k in nzrange(g.A, i))
end

"""
    Graphs.neighbors(g::FactorGraph, f::Factor)

Return a lazy iterator to the edges incident to factor `f`
"""
function Graphs.edges(g::FactorGraph, f::Factor)
    a = f.a
    (FactorGraphEdge(a, g.X.rowval[k], k) for k in nzrange(g.X, a))
end
inedges(g::FactorGraph, x::VariableOrFactor) = Graphs.edges(g, x)
outedges(g::FactorGraph, x::VariableOrFactor) = Graphs.edges(g, x)

"""
    Graphs.neighbors(g::FactorGraph, v::Variable)

Return a lazy iterator to the neighbors of variable `v`
"""
Graphs.neighbors(g::FactorGraph, v::Variable) = @view g.A.rowval[nzrange(g.A, v.i)]

"""
    Graphs.neighbors(g::FactorGraph, f::Factor)

Return a lazy iterator to the neighbors of factor `f`
"""
Graphs.neighbors(g::FactorGraph, f::Factor) = @view g.X.rowval[nzrange(g.X, f.a)]
Graphs.inneighbors(g::FactorGraph, x::VariableOrFactor) = neighbors(g, x)
Graphs.outneighbors(g::FactorGraph, x::VariableOrFactor) = neighbors(g, x)

function Graphs.neighbors(g::FactorGraph, i::Int)
    nvars = nvariables(g)
    if i > nvars
        # i is a factor
        return (j + nvars for j in neighbors(g, Factor(i-nvars)))
    else
        # i is a variable
        return neighbors(g, Variable(i))
    end
end

Graphs.has_vertex(g::FactorGraph, v::Variable) = v.i ≤ nvariables(g)
Graphs.has_vertex(g::FactorGraph, f::Factor) = f.a ≤ nfactors(g)
Graphs.has_vertex(g::FactorGraph, i::Int) = i ≤ nv(g)
Graphs.has_edge(g::FactorGraph, f::Factor, v::Variable) = g.A[f.a, v.i] != 0
Graphs.has_edge(g::FactorGraph, v::Variable, f::Factor) = Graphs.has_edge(g, f, v)

# return true if (src, dst) is an edge in the bipartite view of g
function Graphs.has_edge(g::FactorGraph, src::Integer, dst::Integer)
    Graphs.has_edge(g, Factor(src - nvariables(g)), Variable(dst))
end

Base.zero(g::FactorGraph) = FactorGraph(zero(g.A))

Graphs.is_directed(g::FactorGraph) = false
Graphs.is_directed(::Type{FactorGraph{T}}) where T = false
Graphs.is_bipartite(g::FactorGraph) = true

edge_idx(g::FactorGraph, f::Factor, v::Variable) = nzindex(g.A, f.a, v.i)
edge_idx(g::FactorGraph, v::Variable, f::Factor) = edge_idx(g, f, v)
edge_src_dst(g::FactorGraph, id::Integer) = nzindex(g.A, id)

"""
    get_edge(g::FactorGraph, v::Variable, f::Factor)
    get_edge(g::FactorGraph, f::Factor, v::Variable)
    get_edge(g::FactorGraph, id::Integer)

Get edge given variable and factor or given edge index
""" 
function get_edge(g::FactorGraph, x::VariableOrFactor, y::VariableOrFactor)
    k = edge_idx(g, x, y)
    FactorGraphEdge(x, y, k)
end

function get_edge(g::FactorGraph, id::Integer)
    a, i = edge_src_dst(g, id)
    FactorGraphEdge(a, i, id)
end

"""
    Graphs.LinAlg.adjacency_matrix(g::FactorGraph, T::DataType=Int)

Return the symmetric adjacency matrix of size `nvariables(g) + nfactors(g)` 
where no distinction is made between variable and factor nodes.
Edge indices are preserved in the bottom-left block.
"""
function Graphs.LinAlg.adjacency_matrix(g::FactorGraph, T::DataType=Int)
    m, n = nfactors(g), nvariables(g)
    A = SparseMatrixCSC(g.A.m, g.A.n, g.A.colptr, g.A.rowval, ones(T, nnz(g.A)))
    return [ spzeros(T, n, n)  A'               ;
             A                 spzeros(T, m, m) ]
end

"""
    bipartite_view(g::FactorGraph, T::DataType=Int)

Construct an undirected bipartite `SparseGraph` where vertices are the concatenation
of the variables and factors of the original `FactorGraph`. 
Edge indices are preserved in the bottom-left block.
"""
function bipartite_view(g::FactorGraph, T::DataType=Int)
    A = Graphs.LinAlg.adjacency_matrix(g, T)
    IndexedGraph(A)
end