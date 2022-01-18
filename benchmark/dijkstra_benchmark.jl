using SparseArrays, MatrixNetworks, Graphs, SimpleWeightedGraphs

# using IndexedGraphs
include("../src/IndexedGraphs.jl")
using .IndexedGraphs

# build directed graph
W = sprand(5000, 5000, 0.3)
# remove self loops
for i in 1:size(W,2); W[i,i] = 0; end
dropzeros!(W)
# pick source at random
s = rand(1:size(W,2))

g = IndexedDiGraph(W)
g_MN = MatrixNetwork(W)
g_Graphs = SimpleDiGraph(W)
g_SWG = SimpleWeightedDiGraph(W)

w = nonzeros(permutedims(W))

ds_IG = dijkstra_shortest_paths(g, s, w)
d_IG = ds_IG.dists; p_IG  = ds_IG.parents

d_MN, p_MN = MatrixNetworks.dijkstra(g_MN, s)

ds_Graphs = dijkstra_shortest_paths(g_Graphs, s, W)
d_Graphs = ds_Graphs.dists; p_Graphs = ds_Graphs.parents

ds_SWG = dijkstra_shortest_paths(g_SWG, s, W)
d_SWG = ds_SWG.dists; p_SWG = ds_SWG.parents

# check that results are correct
@assert d_IG == d_MN == d_Graphs == d_SWG
@assert p_IG == p_MN == p_Graphs == p_SWG


### BENCHMARK
using BenchmarkTools
println("IndexedGraphs:")
@btime dijkstra_shortest_paths($g, $s, $w)
println("MatrixNetworks:")
@btime MatrixNetworks.dijkstra($g_MN, $s)
println("Graphs")
@btime dijkstra_shortest_paths($g_Graphs, $s, $W)
println("SimpleWeightedGraphs:")
@btime dijkstra_shortest_paths($g_SWG, $s, $W);