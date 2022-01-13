
using Graphs, SparseArrays

struct SparseMatrixDiGraph{T,TW} <: AbstractGraph{T}
	A::SparseMatrixCSC{Bool, T}
	X::SparseMatrixCSC{T,T}
	W::TW
end

struct IndexedEdge{T}
	src::T
	dst::T
	idx::T
end

function SparseMatrixDiGraph(A::AbstractMatrix{Bool}, W)
	size(A,1) != size(A,2) && throw(ArgumentError("Matrix should be square"))
	any(A[i,i] for i=1:size(A,1)) && throw(ArgumentError("Self loops are not allowed"))
	A = sparse(A)
	length(W) != nnz(A) && throw(ArgumentError("Number of properties different from number of edges"))
	X = sparse(SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, collect(1:nnz(A)))')
	SparseMatrixDiGraph(A, X, W)
end

function SparseMatrixDiGraph(A::AbstractMatrix)
	SparseMatrixDiGraph(SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, fill(true, length(A.nzval))), A.nzval)
end

Graphs.edges(g::SparseMatrixDiGraph) = (IndexedEdge{Int}(i,g.A.rowval[k],k) for i=1:size(g.A,2) for k=nzrange(g.A,i))

Base.eltype(g::SparseMatrixDiGraph{T}) where T = T

Graphs.edgetype(g::SparseMatrixDiGraph{T}) where T = IndexedEdge{T}

Graphs.has_edge(g::SparseMatrixDiGraph, i, j) = g.A[i,j]

Graphs.has_vertex(g::SparseMatrixDiGraph, i) = i ∈ 1:size(g.A, 2)

Graphs.inneighbors(g::SparseMatrixDiGraph, i) = @view g.A.rowval[nzrange(g.A,i)]

Graphs.outneighbors(g::SparseMatrixDiGraph, i) = @view g.X.rowval[nzrange(g.X,i)]

inedges(g::SparseMatrixDiGraph, i) = (IndexedEdge{Int}(i,g.A.rowval[k],k) for k=nzrange(g.A,i))

outedges(g::SparseMatrixDiGraph, i) = (IndexedEdge{Int}(g.X.rowval[k],i,g.X.nzval[k]) for k=nzrange(g.X,i))

property(g::SparseMatrixDiGraph, e::IndexedEdge) = g.W[e.idx]

Graphs.ne(g::SparseMatrixDiGraph) = nnz(g.A)

Graphs.nv(g::SparseMatrixDiGraph) = size(g.A,2)

Graphs.vertices(g::SparseMatrixDiGraph) = 1:size(g.A,2)

Graphs.is_directed(::Type{SparseMatrixDiGraph{T,TW}}) where {T,TW} = true

Graphs.is_directed(g::SparseMatrixDiGraph) = true

Base.zero(g::SparseMatrixDiGraph) = SparseMatrixDiGraph(zero(g.A), similar(g.W, 0))

