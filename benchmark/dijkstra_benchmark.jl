using SparseArrays, MatrixNetworks, Graphs, SimpleWeightedGraphs
using SimpleValueGraphs
# using IndexedGraphs
include("../src/IndexedGraphs.jl")
using .IndexedGraphs
using Random
using LinearAlgebra


# build directed graph
N = 5_000
W = sprand(N, N, 0.3)
# remove self loops
for i in 1:size(W,2); W[i,i] = 0; end
dropzeros!(W)
W = W + W'
w = sparse(UpperTriangular(W)).nzval    # unique wights for Graph
wD = nonzeros(permutedims(W))           # duplicate weights for DiGraph

# pick source at random
rng = MersenneTwister(0)
s = rand(rng, 1:size(W,2))

g = IndexedGraph(W)
gD = IndexedDiGraph(W)
g_MN = MatrixNetwork(W)
g_Graphs = SimpleGraph(W)
g_SWG = SimpleWeightedGraph(W)
g_SVG = ValGraph(g_Graphs, 
    edgeval_types=(Float64,), edgeval_init=(s, d) -> (W[s,d],))


ds_IG = dijkstra_shortest_paths(g, s, w)
d_IG = ds_IG.dists; p_IG  = ds_IG.parents

ds_IDG = dijkstra_shortest_paths(gD, s, wD)
d_IDG = ds_IDG.dists; p_IDG  = ds_IDG.parents

d_MN, p_MN = MatrixNetworks.dijkstra(g_MN, s)

ds_Graphs = dijkstra_shortest_paths(g_Graphs, s, W)
d_Graphs = ds_Graphs.dists; p_Graphs = ds_Graphs.parents

ds_SWG = dijkstra_shortest_paths(g_SWG, s)
d_SWG = ds_SWG.dists; p_SWG = ds_SWG.parents

ds_SVG = SimpleValueGraphs.Experimental.dijkstra_shortest_paths(g_SVG, s)
d_SVG = ds_SVG.dists; p_SVG = ds_SVG.parents

# check that results are correct
@assert d_IG == d_MN == d_Graphs == d_SWG == d_SVG == d_IDG
@assert p_IG == p_MN == p_Graphs == p_SWG == p_SVG == p_IDG


### BENCHMARK
using BenchmarkTools
println("IndexedDiGraph:")
@btime dijkstra_shortest_paths($gD, $s, $wD)
println("IndexedGraph:")
@btime dijkstra_shortest_paths($g, $s, $w)
println("MatrixNetworks:")
@btime MatrixNetworks.dijkstra($g_MN, $s)
println("SimpleGraph")
@btime dijkstra_shortest_paths($g_Graphs, $s, $W)
println("SimpleWeightedGraph:")
@btime dijkstra_shortest_paths($g_SWG, $s);
println("SimpleValueGraph:")
@btime SimpleValueGraphs.Experimental.dijkstra_shortest_paths($g_SVG, $s);