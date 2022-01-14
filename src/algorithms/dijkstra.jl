### Homemade version to compare with MatrixNetworks
using TrackingHeaps

function dijkstra(g::AbstractSparseMatrixGraph{T}, s::Integer, 
        w::AbstractVector{U}=ones(Int, ne(g))) where {T<:Integer, U<:Real}
    n = nv(g)
    inf = typemax(U)
    D = fill(inf, n); D[s] = 0
    parents = fill(0, n); parents[s] = 0
    Q::TrackingHeap{Int64, U, 2, MinHeapOrder, NoTrainingWheels} = TrackingHeaps.TrackingHeap(U; S=TrackingHeaps.NoTrainingWheels)
    sizehint!(Q, n)
    Q[s] = zero(U)
    while !isempty(Q)
        v, d = TrackingHeaps.pop!(Q)
        for e in outedges(g, v)
            i = dst(e)
            dd = d + w[idx(e)]
            if dd < D[i]
                Q[i] = D[i] = dd
                parents[i] = v
            end
        end
    end
    return D, parents
end

##### TEST 
using SparseArrays, MatrixNetworks

# build directed graph
W = sprand(1000, 1000, 0.3)
# remove self loops
for i in 1:size(W,2); W[i,i] = 0; end
dropzeros!(W)
s = rand(1:size(W,2))

g = SparseMatrixDiGraph(W)
g_MN = MatrixNetwork(W)

# w = W.nzval 
w = nonzeros(permutedims(W))

d_SMG, p_SMG = dijkstra(g, s, w)
d_MN, p_MN = MatrixNetworks.dijkstra(g_MN, s)

@assert d_SMG == d_MN
@assert p_SMG == p_MN

using BenchmarkTools
println("SparseMatrixGraphs:")
@btime dijkstra($g, $s, $w)
println("MatrixNetworks:")
@btime MatrixNetworks.dijkstra($g_MN, $s);


### Override for Graphs' version, copied from https://github.com/JuliaGraphs/Graphs.jl/blob/master/src/shortestpaths/dijkstra.jl
using DataStructures, Graphs, SimpleWeightedGraphs
import Graphs.DijkstraState

function Graphs.dijkstra_shortest_paths(g::AbstractSparseMatrixGraph,
    srcs::Vector{U},
    distvec::AbstractVector{T}=ones(Int, ne(g));
    allpaths=false,
    trackvertices=false
    ) where T <: Real where U <: Integer

    nvg = nv(g)
    dists = fill(typemax(T), nvg)
    parents = zeros(U, nvg)

    pathcounts = zeros(nvg)
    preds = fill(Vector{U}(), nvg)
    H = PriorityQueue{U,T}()

    for src in srcs
        dists[src] = zero(T)
        pathcounts[src] = one(Float64)
        H[src] = zero(T)
    end

    closest_vertices = Vector{U}()  # Maintains vertices in order of distances from source
    sizehint!(closest_vertices, nvg)

    while !isempty(H)
        u = dequeue!(H)

        if trackvertices
            push!(closest_vertices, u)
        end

        d = dists[u] # Cannot be typemax if `u` is in the queue
        for e in outedges(g, u)
            alt = d + distvec[idx(e)]
            v = dst(e)

            visited = dists[v] != typemax(T)
            if !visited
                dists[v] = alt
                parents[v] = u

                pathcounts[v] += pathcounts[u]
                if allpaths
                    preds[v] = [u;]
                end
                H[v] = alt
            elseif alt < dists[v]
                dists[v] = alt
                parents[v] = u
                #615
                pathcounts[v] = pathcounts[u]
                if allpaths
                    resize!(preds[v], 1)
                    preds[v][1] = u
                end
                H[v] = alt
            elseif alt == dists[v]
                pathcounts[v] += pathcounts[u]
                if allpaths
                    push!(preds[v], u)
                end
            end
        end
    end

    if trackvertices
        for s in vertices(g)
            visited = dists[s] != typemax(T)
            if !visited[s]
                push!(closest_vertices, s)
            end
        end
    end

    for src in srcs
        pathcounts[src] = one(Float64)
        parents[src] = 0
        empty!(preds[src])
    end

    return DijkstraState{T,U}(parents, dists, preds, pathcounts, closest_vertices)
end

Graphs.dijkstra_shortest_paths(g::AbstractSparseMatrixGraph, src::Integer, distvec::AbstractVector=ones(Int, ne(g)); allpaths=false, trackvertices=false) =
Graphs.dijkstra_shortest_paths(g, [src;], distvec; allpaths=allpaths, trackvertices=trackvertices)

g_Graphs = SimpleWeightedDiGraph(W)

ds_SMG = Graphs.dijkstra_shortest_paths(g, s, w)
ds_Graphs = Graphs.dijkstra_shortest_paths(g_Graphs, s)

@assert ds_SMG.dists == ds_Graphs.dists
@assert ds_SMG.parents == ds_Graphs.parents

println()
println("SparseMatrixGraphs:")
@btime Graphs.dijkstra_shortest_paths($g, $s, $w)
println("SimpleWeightedGraphs:")
@btime Graphs.dijkstra_shortest_paths($g_Graphs, $s);