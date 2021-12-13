using Graphs, SparseArrays, LinearAlgebra

@testset "graphs" begin
	A = sparse(I(10))
	B = spdiagm(1=>trues(9),-1=>trues(9))
	gA = SparseMatrixDiGraph(A)
	gB = SparseMatrixDiGraph(B)
	@test is_directed(gA)
	@test !is_connected(gA)
	@test is_connected(gB)
end
