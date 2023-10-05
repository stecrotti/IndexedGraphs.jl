using Graphs, IndexedGraphs, SparseArrays

@testset "bipartite graph" begin
    nl = 15
    nr = 27
    n = nl + nr
    g = complete_bipartite_graph(nl, nr) |> IndexedGraph
    gb = BipartiteIndexedGraph(g)
    @test vertices(g) == vertices(gb)
    @test collect(edges(g)) == collect(edges(gb))
    @test nv(g) == nv(gb)
    @test ne(g) == ne(gb)
    @test all(collect(neighbors(g,i)) == collect(neighbors(gb, i)) for i in vertices(g))
    @test all(collect(inedges(g,i)) == collect(inedges(gb, i)) for i in vertices(g))
    @test all(collect(outedges(g,i)) == collect(outedges(gb, i)) for i in vertices(g))

    distmx = adjacency_matrix(g) .* rand(n,n)
    distvec = nonzeros(permutedims(distmx))
    sources = rand(1:n, 10)
    d = dijkstra_shortest_paths(g, sources, distvec)
    db = dijkstra_shortest_paths(gb, sources, distvec)
    @test all(getproperty(d, p) == getproperty(db, p) for p in fieldnames(typeof(d)))
end