# IndexedGraphs.jl
Not all graph edges come with an index. These do.

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
