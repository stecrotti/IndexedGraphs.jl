using SparseArrays, Graphs

A = sprand(Bool, 20, 20, 0.5)
for i in 1:20; A[i,i] = 0; end
dropzeros!(A)
g = IndexedDiGraph(A)

@testset "directed graph" begin

    @testset "transpose constructor" begin
        At = sparse(A')
        gg = IndexedDiGraph( transpose(At) )
        @test gg.A.rowval === At.rowval
        @test gg.A.rowval == g.A.rowval
    end

    @testset "basics" begin
        @test is_directed(g)
        @test length(collect(edges(g))) == ne(g)
        i = 3
        oute = outedges(g, i)
        outn = outneighbors(g, i)
        @test all(dst(e) == j for (e,j) in zip(oute, outn))
        @test all(src(e) == i for e in oute)
    end

    @testset "edge indexing" begin
        @test all( e == get_edge(g, src(e), dst(e)) for e in edges(g) )
        @test all( e == get_edge(g, idx(e)) for e in edges(g) )

        passed = falses(ne(g))
        for (i,e) in enumerate(edges(g))
            id = idx(get_edge(g, src(e), dst(e)))
            ee = get_edge(g, id)
            passed[i] = ee == e
        end
        @test all(passed)
    end

    @testset "construct from SimpleGraph" begin
        sg = SimpleDiGraph(A)
        ig = IndexedDiGraph(sg)
        @test adjacency_matrix(sg) == adjacency_matrix(ig)
        S = A + A'
        sg = SimpleDiGraph(S)
        ig = IndexedDiGraph(sg)
        @test adjacency_matrix(sg) == adjacency_matrix(ig)
    end
end