using SparseArrays, Graphs, IndexedGraphs, LinearAlgebra

function dotest(g, g_Graphs, s, w, W, allpaths, trackvertices)
    ds_IG = dijkstra_shortest_paths(g, s, w; allpaths, trackvertices)
    ds_Graphs = dijkstra_shortest_paths(g_Graphs, s, W; allpaths, trackvertices)
    # check that results are correct
    @test ds_IG.dists == ds_Graphs.dists
    @test ds_IG.parents == ds_Graphs.parents
    @test ds_IG.predecessors == ds_Graphs.predecessors
    @test ds_IG.pathcounts == ds_Graphs.pathcounts
    @test ds_IG.closest_vertices == ds_Graphs.closest_vertices
end

@testset "dijkstra " begin
    @testset "directed" begin
        for G in (IndexedDiGraph, IndexedBiDiGraph)
            for allpaths = (true,false)
                for trackvertices = (true, false)
                    for N=(10,20,30)
                        # build directed graph
                        W = sprand(N, N, 0.5)
                        # remove self loops
                        for i in 1:size(W,2); W[i,i] = 0; end
                        dropzeros!(W)
                        # pick sources at random
                        s = rand(1:size(W,2), 2)

                        g = G(W)
                        g_Graphs = SimpleDiGraph(W)

                        w = nonzeros(permutedims(W))
                        dotest(g, g_Graphs, s, w, W, allpaths, trackvertices)
                    end
                end
            end
        end
    end
end
