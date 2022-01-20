# IndexedBiDiGraph

A type representing directed graphs. 

Use this when you need to access both inedges and outedges (or inneighbors and outneighbors).
For a lighter data structure check out [IndexedDiGraph](@ref).

```@docs
IndexedBiDiGraph
```

```@docs
IndexedBiDiGraph(A::AbstractMatrix)
```
Example:
```@example
using SparseArrays, IndexedGraphs
At = sprand(100, 100, 0.1)           # At[i,j] corresponds to edge j=>i
for i in 1:100; At[i,i] = 0; end
dropzeros!(At)
g = IndexedBiDiGraph(transpose(At))  
g.A.rowval === At.rowval
```


