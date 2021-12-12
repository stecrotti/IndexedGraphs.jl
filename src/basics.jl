
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
	A = sparse(A)
	X = sparse(SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, collect(1:nnz(A)))')
	SparseMatrixDiGraph(A, X, W)
end

function SparseMatrixDiGraph(A::AbstractMatrix)
	SparseMatrixDiGraph(SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, fill(true, length(A.nzval))), A.nzval)
end

edges(g::SparseMatrixDiGraph) = (IndexedEdge{Int}(i,g.A.rowval[k],k) for i=1:size(g.A,2) for k=nzrange(g.A,i))

Base.eltype(g::SparseMatrixDiGraph{T}) where T= T

edgetype(g::SparseMatrixDiGraph{T}) where T = IndexedEdge{T}

has_edge(g::SparseMatrixDiGraph, i, j) = g.A[i,j]

has_vertex(g::SparseMatrixDiGraph, i) = 1 ≤ i ≤ size(g.A, 2)

inneighbors(g::SparseMatrixDiGraph, i) = @view g.A.rowval[nzrange(g.A,i)]

outneighbors(g::SparseMatrixDiGraph, i) = @view g.A.rowval[g.X.nzval[nzrange(g.X,i)]]

outedges(g::SparseMatrixDiGraph, i) = (IndexedEdge{Int}(i,g.A.rowval[k],k) for k=nzrange(g.A,i))

inedges(g::SparseMatrixDiGraph, i) = (IndexedEdge{Int}(g.X.rowval[k],i,k) for k=@view g.X.nzval[nzrange(g.X,i)])

property(g::SparseMatrixDiGraph, e::IndexedEdge) = g.W[e.idx]

ne(g::SparseMatrixDiGraph) = nnz(g.A)

nv(g::SparseMatrixDiGraph) = size(A,2)

vertices(g::SparseMatrixDiGraph) = 1:size(A,2)

is_directed(g::SparseMatrixDiGraph) = true

zero(g::SparseMatrixDiGraph) = SparseMatrixDiGraph(zero(A), @view W[1:0])
