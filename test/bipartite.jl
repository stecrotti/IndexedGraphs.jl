using Graphs, IndexedGraphs, SparseArrays, Random

@testset "bipartite graph" begin
    nl = 15
    nr = 27
    n = nl + nr
    g = complete_bipartite_graph(nl, nr) |> IndexedGraph
    gb = BipartiteIndexedGraph(g)

    @testset "static properties" begin
        @test !is_directed(gb)
        @test !is_directed(typeof(gb))
        @test is_bipartite(gb)
    end

    @testset "basics" begin
        @test vertices(g) == vertices(gb)
        @test collect(edges(g)) == collect(edges(gb))
        @test nv(g) == nv(gb)
        @test ne(g) == ne(gb)
        @test all(collect(neighbors(g,i)) == collect(neighbors(gb, i)) for i in vertices(g))
        @test all(collect(inedges(g,i)) == collect(inedges(gb, i)) for i in vertices(g))
        @test all(collect(outedges(g,i)) == collect(outedges(gb, i)) for i in vertices(g))
        @test all(has_edge(gb, s, d) for (s, d) in edges(gb))
        @test all(degree(g, i) == length(collect(neighbors(g, i))) for i in vertices(g))
        @test all(degree(g, i) == length(collect(inedges(g, i))) for i in vertices(g))
        @test adjacency_matrix(gb) == adjacency_matrix(g)
    end

    @testset "degree" begin
        @test all(1:nv_left(gb)) do iL
            v = vertex(iL, Left)
            i = linearindex(gb, v)
            degree(gb, v) == degree(g, i)
        end
        @test all(1:nv_right(gb)) do iR
            v = vertex(iR, Right)
            i = linearindex(gb, v)
            degree(gb, v) == degree(g, i)
        end
    end

    @testset "left and right" begin
        vl = [linearindex(gb, i, Left) for i in 1:nv_left(gb)]
        vr = [linearindex(gb, i, Right) for i in 1:nv_right(gb)]
        @test vl == vertices_left(gb)
        @test vr == vertices_right(gb)

        @test all(all(vertex_left(gb, e)==l for e in outedges(gb, vertex(l, Left))) 
            for l in 1:nv_left(gb))
        @test all(all(vertex_left(gb, e)==l for e in inedges(gb, vertex(l, Left))) 
            for l in 1:nv_left(gb))
        @test all(all(vertex_right(gb, e)==r for e in outedges(gb, vertex(r, Right))) 
            for r in 1:nv_right(gb))
        @test all(all(vertex_right(gb, e)==r for e in inedges(gb, vertex(r, Right))) 
            for r in 1:nv_right(gb))
    end

    @testset "dijkstra" begin
        distmx = adjacency_matrix(g) .* rand(n,n)
        distvec = nonzeros(permutedims(distmx))
        sources = rand(1:n, 10)
        d = dijkstra_shortest_paths(g, sources, distvec)
        db = dijkstra_shortest_paths(gb, sources, distvec)
        @test all(getproperty(d, p) == getproperty(db, p) for p in fieldnames(typeof(d)))
    end

    @testset "bipartite generators" begin
        rng = MersenneTwister(0)
        ngraphs = 20
        nrights = rand(rng, 5:50, ngraphs)
        nlefts = rand(rng, 5:50, ngraphs)
        es = [rand(rng, 1:n*m) for (n, m) in zip(nrights, nlefts)]
        
        @testset "Random bipartite graph - fixed # edges" begin
            @test all(zip(nrights, nlefts, es)) do (n, m, e)
                g = rand_bipartite_graph(rng, m, n, e)
                nv_right(g) == n && nv_left(g) == m && ne(g) == e
            end
        end
        
        @testset "Random bipartite graph - prob of edges" begin
            p = 0.1
            @test all(zip(nrights, nlefts)) do (n, m)
                g = rand_bipartite_graph(rng, n, m, p)
                nv_right(g) == n && nv_left(g) == m
            end
        end
        
        @testset "Random regular bipartite graph" begin
            k = 4
            @test all(zip(nrights, nlefts)) do (n, m)
                g = rand_regular_bipartite_graph(rng, n, m, k)
                nv_right(g) == n && nv_left(g) == m && ne(g) == m * k
            end
        end
        
        @testset "Random bipartite tree" begin
            @test all(nrights) do n
                g = rand_bipartite_tree(rng, n)
                nv(g) == n && !is_cyclic(g)
            end
        end
    end
end