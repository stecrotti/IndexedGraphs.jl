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

# W = spzeros(10, 10)
# W[2,3] = W[3,4] = W[4,5] = 1
# s = 2

g = SparseMatrixDiGraph(W)
g_MN = MatrixNetwork(W)

# w = W.nzval 
w = nonzeros(permutedims(W))

d_SMG, p_SMG = dijkstra(g, s, w)
d_MN, p_MN = MatrixNetworks.dijkstra(g_MN, s)

@assert d_SMG == d_MN
@assert p_SMG == p_MN

using BenchmarkTools
@btime dijkstra($g, $s, $w)
@btime MatrixNetworks.dijkstra($g_MN, $s);