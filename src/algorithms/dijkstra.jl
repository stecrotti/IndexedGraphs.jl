### Override for Graphs' version, copied from https://github.com/JuliaGraphs/Graphs.jl/blob/master/src/shortestpaths/dijkstra.jl with some adjustments
function dijkstra_shortest_paths(g::AbstractIndexedGraph,
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
    H::TrackingHeap{Int64, T, 2, MinHeapOrder, NoTrainingWheels} = TrackingHeap(T; S=NoTrainingWheels)

    for src in srcs
        dists[src] = zero(T)
        pathcounts[src] = one(Float64)
        H[src] = zero(T)
    end

    closest_vertices = Vector{U}()  # Maintains vertices in order of distances from source
    sizehint!(closest_vertices, nvg)

    while !isempty(H)
        u, d = pop!(H)

        if trackvertices
            push!(closest_vertices, u)
        end

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

function dijkstra_shortest_paths(g::AbstractIndexedGraph, src::Integer, 
        distvec::AbstractVector=ones(Int, ne(g)); kw...)
    dijkstra_shortest_paths(g, [src;], distvec; kw...)
end

function dijkstra_shortest_paths(g::AbstractIndexedGraph, src, 
        distmx::AbstractMatrix; kw...)
    Wt = convert(SparseMatrixCSC, transpose(distmx))
    distvec = nonzeros(Wt)
    dijkstra_shortest_paths(g, src, distvec; kw...)
end


### Homemade minimal version
function dijkstra(g::AbstractIndexedGraph{T}, s::Integer, 
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