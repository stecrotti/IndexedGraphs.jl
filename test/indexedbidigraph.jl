using SparseArrays, Graphs

A = sprand(Bool, 20, 20, 0.5)
for i in 1:20; A[i,i] = 0; end
dropzeros!(A)
g = IndexedBiDiGraph(A)

@testset "BiDirected graph" begin

    @testset "show" begin
        buf = IOBuffer()
        show(buf, g)
        @test String(take!(buf)) == "{20, $(ne(g))} IndexedBiDiGraph{$(Int)}\n"
    end

    @testset "transpose constructor" begin
        At = sparse(A')
        gg = IndexedBiDiGraph( transpose(At) )
        @test gg.A.rowval === At.rowval
        @test gg.A.rowval == g.A.rowval
    end

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
        ig = IndexedBiDiGraph(sg)
        @test adjacency_matrix(sg) == adjacency_matrix(ig)
        S = A + A'
        sg = SimpleDiGraph(S)
        ig = IndexedBiDiGraph(sg)
        @test adjacency_matrix(sg) == adjacency_matrix(ig)
    end
end