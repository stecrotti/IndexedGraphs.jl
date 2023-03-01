using SparseArrays, Graphs

A = sprand(Bool, 20, 20, 0.5)
for i in 1:20; A[i,i] = 0; end
A = A + A'
dropzeros!(A)
g = IndexedGraph(A)

@testset "undirected graph" begin

    @testset "show" begin
        buf = IOBuffer()
        show(buf, g)
        @test String(take!(buf)) == "{20, $(ne(g))} undirected IndexedGraph{$(Int)}\n"
    end

    @testset "basics" begin
        @test !is_directed(g)
        @test length(collect(edges(g))) == ne(g)
        i = 3
        es = inedges(g, i)
        neigs = neighbors(g, i)
        @test all(
            ((src(e), dst(e)) == (i, j)) || ((src(e), dst(e)) == (j, i))
            for (e, j) in zip(es, neigs)
        )
    end

    @testset "edge indexing" begin
        @test all( e == get_edge(g, src(e), dst(e)) for e in edges(g) )
        @test all( e == get_edge(g, idx(e)) for e in edges(g) )
    end

    @testset "construct from SimpleGraph" begin
        sg = SimpleGraph(A)
        ig = IndexedGraph(sg)
        @test adjacency_matrix(sg) == adjacency_matrix(ig)
    end
end

@testset "iterator within edge" begin
    es = map(edges(g)) do e
        i, j, ij = e
        i==src(e) && j == dst(e) && ij == idx(e)
    end
    @test all(es)    
end