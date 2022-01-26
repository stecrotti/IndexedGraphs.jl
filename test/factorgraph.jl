using SparseArrays, Graphs

A = sprand(Bool, 10, 20, 0.5)
for i in 1:10; A[i,i] = 0; end
dropzeros!(A)
g = FactorGraph(A)

@testset "factor graph" begin
    @testset "basics" begin
        @test is_bipartite(g)
        @test !is_directed(g)
        @test length(collect(edges(g))) == ne(g)
        @test length(collect(vertices(g))) == nv(g)
    end

    @testset "neighbors of variable" begin
        i = 3
        es = edges(g, Variable(i))
        neigs = neighbors(g, Variable(i))
        @test all( (e.a, e.i) == (a, i) for (e, a) in zip(es, neigs) )
    end

    @testset "neighbors of factor" begin
        a = 9
        es = edges(g, Factor(a))
        neigs = neighbors(g, Factor(a))
        @test all( (e.a, e.i) == (a, i) for (e, i) in zip(es, neigs) )
    end

    @testset "edge indexing" begin
        @test all( e == get_edge(g, Factor(e.a), Variable(e.i)) for e in edges(g) )
        @test all( e == get_edge(g, idx(e)) for e in edges(g) )
    end

    @testset "bipartiteness" begin
        gg = bipartite_view(g)
        @test is_bipartite(gg)
    end
end