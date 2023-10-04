using SparseArrays, Graphs

A = sprand(Bool, 10, 20, 0.5)
for i in 1:10; A[i,i] = 0; end
dropzeros!(A)
g = FactorGraph(A)

@testset "factor graph" begin

    @testset "show" begin
        buf = IOBuffer()
        show(buf, g)
        @test String(take!(buf)) == "FactorGraph{$(Int)} with $(nvariables(g)) variables,"*
            " $(nfactors(g)) factors\n"
    end

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

    @testset "from pairwise graph" begin
        B = sprand(Bool, 20, 20, 0.5)
        for i in 1:20; B[i,i] = 0; end
        B = B + B'
        dropzeros!(B)
        gpair = IndexedGraph(B)
        gfact = FactorGraph(gpair)
        @test nfactors(gfact) == ne(gpair)
        @test nvariables(gfact) == nv(gpair)
        @test all(
            has_edge(gfact, Variable(i), Factor(id)) &&
            has_edge(gfact, Variable(j), Factor(id))
            for (i, j, id) in edges(gpair)
            )
    end
end