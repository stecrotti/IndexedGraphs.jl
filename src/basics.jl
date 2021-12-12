
using Graphs, SparseArrays

struct MatrixDiGraph{T,TW} <: AbstractGraph{T}
	A::SparseMatrixCSC{Bool, T}
	X::SparseMatrixCSC{T,T}
	W::TW
end

struct IndexedEdge{T}
	src::T
	dst::T
	idx::T
end

function MatrixDiGraph(A::AbstractMatrix{Bool}, W)
	A = sparse(A)
	X = sparse(SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, collect(1:nnz(A)))')
	MatrixDiGraph(A, X, W)
end

function MatrixDiGraph(A::AbstractMatrix)
	MatrixDiGraph(SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, fill(true, length(A.nzval))), A.nzval)
end

edges(g::MatrixDiGraph) = (IndexedEdge{Int}(i,g.A.rowval[k],k) for i=1:size(g.A,2) for k=nzrange(g.A,i))

Base.eltype(g::MatrixDiGraph{T}) where T= T

edgetype(g::MatrixDiGraph{T}) where T = IndexedEdge{T}

has_edge(g::MatrixDiGraph, i, j) = g.A[i,j]

has_vertex(g::MatrixDiGraph, i) = 1 ≤ i ≤ size(g.A, 2)

inneighbors(g::MatrixDiGraph, i) = @view g.A.rowval[nzrange(g.A,i)]

outneighbors(g::MatrixDiGraph, i) = @view g.A.rowval[g.X.nzval[nzrange(g.X,i)]]

outedges(g::MatrixDiGraph, i) = (IndexedEdge{Int}(i,g.A.rowval[k],k) for k=nzrange(g.A,i))

inedges(g::MatrixDiGraph, i) = (IndexedEdge{Int}(g.X.rowval[k],i,k) for k=@view g.X.nzval[nzrange(g.X,i)])

weight(g::MatrixDiGraph, e::IndexedEdge) = g.W[e.idx]

ne(g::MatrixDiGraph) = nnz(g.A)

nv(g::MatrixDiGraph) = size(A,2)

vertices(g::MatrixDiGraph) = 1:size(A,2)

is_directed(g::MatrixDiGraph) = true

zero(g::MatrixDiGraph) = MatrixDiGraph(zero(A), @view W[1:0])
