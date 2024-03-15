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

@testset "SymmetricBiDirected graph" begin

    A = sprand(Bool, 20, 20, 0.5)
    for i in 1:20; A[i,i] = 0; end
    A = A + A'
    dropzeros!(A)
    g = IndexedGraph(A)
    gd, dir2undir, undir2dir = bidirected_with_mappings(g)
    gbd = IndexedBiDiGraph(A)

    @testset "mappings" begin
        eu = edges(g) |> collect    # undirected edges
        ed = edges(gd) |> collect   # directed edges
        @test all( let
            e = eu[dir2undir[idd]]
            src(e) == min(i,j) && dst(e) == max(i,j)
        end for (i,j,idd) in ed
        )
        @test all( let
            es = ed[undir2dir[idu]]
            src(es[1]) == dst(es[2]) && src(es[2]) == dst(es[1])
        end for (i,j,idu) in eu
        )
    end

    @testset "show" begin
        buf = IOBuffer()
        show(buf, gd)
        @test String(take!(buf)) == "{20, $(ne(gd))} SymmetricIndexedBiDiGraph{$(Int)}\n"
    end

    @testset "basics" begin
        @test is_directed(typeof(gbd)) == is_directed(typeof(gd))
        @test is_directed(gbd) == is_directed(gd)
        @test length(collect(edges(gbd))) == ne(gbd) == length(collect(edges(gd))) == ne(gd)
        i = 3
        ine = inedges(gbd, i); ined = inedges(gd, i)
        inn = inneighbors(gbd, i); innd = inneighbors(gd, i)
        @test all(src(e) == j == src(ed) == jd for (e,j,ed,jd) in zip(ine, inn, ined, innd))
        @test all(dst(e) == i == dst(ed) for (e, ed) in zip(ine, ined))
    end

    @testset "edge indexing" begin
        @test all( e == get_edge(gbd, src(e), dst(e)) for e in edges(gd) )
        @test all( e == get_edge(gbd, idx(e)) for e in edges(gd) )
        @test all(let
            id = idx(get_edge(gd, src(e), dst(e)))
            ee = get_edge(gd, id)
            ee == e
        end for (i,e) in enumerate(edges(gd)))
    end

    @testset "Graphs.jl methods" begin
        @test edgetype(gd) == edgetype(gbd)
        @test is_bipartite(gd) == is_bipartite(gbd)
        @test adjacency_matrix(gd) == adjacency_matrix(gbd)
        @test all(let
            has_vertex(gd, i) == has_vertex(gbd, i) &&
            neighbors(gd, i) == neighbors(gbd, i) &&
            outneighbors(gd, i) == outneighbors(gbd, i) &&
            degree(gd, i) == degree(gbd, i) &&
            collect(outedges(gd, i)) == collect(outedges(gbd, i))
        end for i in 1:20)
        @test all(has_edge(gd, i, j) == has_edge(gbd, i, j) for i in 1:20 for j in 1:20)
    end
end