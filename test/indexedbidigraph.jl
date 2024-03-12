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
        @test is_directed(typeof(g))
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

@testset "complete BiDirected graph" begin

    A = sprand(Bool, 20, 20, 0.5)
    for i in 1:20; A[i,i] = 0; end
    A = A + A'
    dropzeros!(A)
    g = CompleteIndexedBiDiGraph(IndexedGraph(A))
    gd = IndexedBiDiGraph(A)

    @testset "show" begin
        buf = IOBuffer()
        show(buf, g)
        @test String(take!(buf)) == "{20, $(ne(g))} CompleteIndexedBiDiGraph{$(Int)}\n"
    end

    @testset "basics" begin
        @test is_directed(typeof(g)) == is_directed(typeof(gd))
        @test is_directed(g) == is_directed(gd)
        @test length(collect(edges(g))) == ne(g) == length(collect(edges(gd))) == ne(gd)
        i = 3
        ine = inedges(g, i); ined = inedges(gd, i)
        inn = inneighbors(g, i); innd = inneighbors(gd, i)
        @test all(src(e) == j == src(ed) == jd for (e,j,ed,jd) in zip(ine, inn, ined, innd))
        @test all(dst(e) == i == dst(ed) for (e, ed) in zip(ine, ined))
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
end