# IndexedGraphs.jl

_Not all edges come with an index. These do_

[![CI](https://github.com/stecrotti/IndexedGraphs.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/stecrotti/IndexedGraphs.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/stecrotti/IndexedGraphs.jl/branch/main/graph/badge.svg?token=CYLRPHU098)](https://codecov.io/gh/stecrotti/IndexedGraphs.jl)

## Documentation
https://stecrotti.github.io/IndexedGraphs.jl/dev

## Overview
A **Graphs.jl**-compatible implementation of [SparseMatrixCSC](https://github.com/JuliaLang/SparseArrays.jl)-based graphs, allowing fast access to arbitrary edge properties

* The code implements the **Graphs.jl** interface for directed and undirected graphs.
* Edge properties live separate from the graph, so different sets of properties can be associated to the same graph.
* In addition, it implements `inedges` and `outedges` for O(1) access to neighborhood
* Edges are indexed, and the index can be used to access edge properties very efficiently.
* `IndexedBiDirectedGraphs` store both the direct and the transposed adjancency matrix for efficient access

A number of other packages implement graph based on CSC matrix representation or similar, namely **StaticGraphs**, **SimpleWeightedGraphs** and **MatrixNetworks**

* [StaticGraphs](https://github.com/JuliaGraphs/StaticGraphs.jl): No edge properties
* [SimpleWeightedGraphs](https://github.com/JuliaGraphs/SimpleWeightedGraphs.jl): Also based on `SparseMatrixCSC`, allows for numerical edge properties. However, no edge weight can be 0 (or otherwise the edge is sometimes removed), and does not allow arbitrary edge properties
* [MatrixNetworks](https://github.com/JuliaGraphs/MatrixNetworks.jl): Also based on `SparseMatrixCSC`, allows for numerical edge properties. However, no edge weight can be 0 (or otherwise the edge is sometimes removed), and does not allow arbitrary edge properties. Does not implement the Graphs interface.

## Navigating graphs
The most natural and efficient way to iterate over an `IndexedGraph` is to iterate over neighboring nodes or edges
```
julia> A = [0 0 1;
            1 0 0;
            1 1 0];
julia> g = IndexedDiGraph(A);
julia> i = 3;
julia> out_i = outedges(g, i);
julia> collect(out_i)
2-element Vector{IndexedGraphs.IndexedEdge{Int64}}:
 Indexed Edge 3 => 1 with index 4
 Indexed Edge 3 => 2 with index 5
```
Edge indices, 4 and 5 in this case, can be extracted with `idx` and used to access properties stored in a separate container
```
julia> e = first(out_i);
julia> src(e), dst(e), idx(e)
(3, 1, 4)
```

## Benchmark
Performance on [Dijkstra's algorithm](https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm) compared with the packages listed above, as computed [here](https://github.com/stecrotti/IndexedGraphs.jl/blob/main/benchmark/dijkstra_benchmark.jl) for a random directed graph of size 5000.
```
IndexedGraphs:
  28.865 ms (36 allocations: 535.31 KiB)
MatrixNetworks:
  29.275 ms (14 allocations: 1.09 MiB)
Graphs
  2.034 s (45 allocations: 808.39 KiB)
SimpleWeightedGraphs:
  1.903 s (45 allocations: 808.39 KiB)
```
