using SparseArrays, Graphs

A = sprand(Bool, 20, 20, 0.5)
for i in 1:20; A[i,i] = 0; end
dropzeros!(A)
g = SparseMatrixDiGraph(A)

@testset "directed graph" begin
    @testset "basics" begin
        @test is_directed(g)
        @test length(collect(edges(g))) == ne(g)
        i = 3
        ine = inedges(g, i)
        inn = inneighbors(g, i)
        @test all(src(e) == j for (e,j) in zip(ine, inn))
        @test all(dst(e) == i for e in ine)
    end

    @testset "edge indexing" begin
        for e in edges(g)
            @test e == get_edge(g, src(e), dst(e)) 
            @test e == get_edge(g, idx(e))
        end
    end
end