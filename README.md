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
```julia
A = [0 0 1;
     1 0 0;
     1 1 0]
g = IndexedDiGraph(A)
i = 3
out_i = outedges(g, i)
collect(out_i)
```
outputs:
```julia
2-element Vector{IndexedGraphs.IndexedEdge{Int64}}:
 Indexed Edge 3 => 1 with index 4
 Indexed Edge 3 => 2 with index 5
```
Edge indices, 4 and 5 in this case, can be extracted with `idx` and used to access properties stored in a separate container
```julia
e = first(out_i)
src(e), dst(e), idx(e)
```
outputs:
```julia
(3, 1, 4)
```

## Benchmark
Performance on [Dijkstra's algorithm](https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm) compared with the packages listed above, as computed [here](https://github.com/stecrotti/IndexedGraphs.jl/blob/main/benchmark/dijkstra_benchmark.jl) for a random symmetric weight matrix with 10^4 nodes and ~10^5 edges.

```julia
IndexedDiGraph:
  2.840 ms (22 allocations: 547.91 KiB)
IndexedGraph:
  3.131 ms (22 allocations: 547.91 KiB)
MatrixNetwork:
  3.031 ms (13 allocations: 407.45 KiB)
SimpleGraph
  11.935 ms (45 allocations: 1008.58 KiB)
SimpleWeightedGraph:
  10.610 ms (45 allocations: 1008.58 KiB)
ValGraph (SimpleValueGraphs.Experimental):
  6.620 ms (48 allocations: 1000.06 KiB)
```

**Note**: For an undirected graph, `IndexedGraph` gives one unique index to each undirected edge (`i=>j` and `j=>i` have the same index). This makes the memory layout less efficient when traversing the graph (although it is very efficient to modify the properties compared with the alternatives). 
If no property modification is needed, as is the case with Dijkstra, it is more convenient to just employ an `IndexedDiGraph` with symmetric edges and weights.
