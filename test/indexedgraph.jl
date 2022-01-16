using SparseArrays, Graphs

A = sprand(Bool, 20, 20, 0.5)
for i in 1:20; A[i,i] = 0; end
A = A + A'
dropzeros!(A)
g = IndexedGraph(A)

@testset "undirected graph" begin
    @testset "basics" begin
        @test !is_directed(g)
        @test length(collect(edges(g))) == ne(g)
        i = 3
        es = inedges(g, i)
        neigs = neighbors(g, i)
        for (e, j) in zip(es, neigs)
            @test ((src(e), dst(e)) == (i, j)) || ((src(e), dst(e)) == (j, i))
        end
    end

    @testset "edge indexing" begin
        for e in edges(g)
            @test e == get_edge(g, src(e), dst(e)) 
            @test e == get_edge(g, idx(e))
        end
    end
end