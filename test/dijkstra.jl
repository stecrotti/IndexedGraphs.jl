using SparseArrays, Graphs, IndexedGraphs, LinearAlgebra

function dotest(ds_IG, ds_Graphs)
    # check that results are correct
    @test ds_IG.dists == ds_Graphs.dists
    @test ds_IG.parents == ds_Graphs.parents
    @test ds_IG.predecessors == ds_Graphs.predecessors
    @test ds_IG.pathcounts == ds_Graphs.pathcounts
    @test ds_IG.closest_vertices == ds_Graphs.closest_vertices
end

@testset "dijkstra " begin
    @testset "directed" begin
        for allpaths = (true,false)
            for trackvertices = (true, false)
                for N=(10,20,30)
                    # build directed graph
                    W = sprand(N, N, 0.5)
                    # remove self loops
                    for i in 1:size(W,2); W[i,i] = 0; end
                    dropzeros!(W)
                    # pick sources at random
                    s = rand(1:N, 2)
                    g_Graphs = SimpleDiGraph(W)
                    ds_Graphs = dijkstra_shortest_paths(g_Graphs, s, W; allpaths=allpaths, trackvertices=trackvertices)
                    for G in (IndexedDiGraph, IndexedBiDiGraph)
                        g = G(W)
                        w = nonzeros(permutedims(W))
                        ds_IG = dijkstra_shortest_paths(g, s, w; allpaths=allpaths, trackvertices=trackvertices)
                        dotest(ds_IG, ds_Graphs)
                    end
                end
            end
        end
    end
    @testset "undirected" begin
        for allpaths = (true,false)
            for trackvertices = (true, )
                for N = (10,20,30)
                    W = sprand(N, N, 0.5)
                    W[diagind(W)] .= 0
                    Wl = triu(W)
                    dropzeros!(Wl)
                    W = Wl .+ Wl'
                    g = IndexedGraph(W)
                    g_Graphs = SimpleDiGraph(W)
                    s = rand(1:N, 2)
                    ds_Graphs = dijkstra_shortest_paths(g_Graphs, s, W; allpaths=allpaths, trackvertices=trackvertices)
                    ds_IG = dijkstra_shortest_paths(g, s, Wl.nzval; allpaths=allpaths, trackvertices=trackvertices)
                    dotest(ds_IG, ds_Graphs)
                end
            end
        end
    end
end
