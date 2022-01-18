var documenterSearchIndex = {"docs":
[{"location":"graph/#IndexedGraph","page":"IndexedGraph","title":"IndexedGraph","text":"","category":"section"},{"location":"graph/","page":"IndexedGraph","title":"IndexedGraph","text":"A type representing undirected graphs.","category":"page"},{"location":"graph/","page":"IndexedGraph","title":"IndexedGraph","text":"IndexedGraph","category":"page"},{"location":"graph/#IndexedGraphs.IndexedGraph","page":"IndexedGraph","title":"IndexedGraphs.IndexedGraph","text":"IndexedGraph{T<:Integer} <: AbstractIndexedGraph{T}\n\nA type representing a sparse undirected graph.\n\nFIELDS\n\nA – square adjacency matrix. A[i,j] == A[j,i] contains the unique index associated to unidrected edge (i,j)\n\n\n\n\n\n","category":"type"},{"location":"graph/","page":"IndexedGraph","title":"IndexedGraph","text":"IndexedGraph(A::AbstractMatrix)","category":"page"},{"location":"graph/#IndexedGraphs.IndexedGraph-Tuple{AbstractMatrix}","page":"IndexedGraph","title":"IndexedGraphs.IndexedGraph","text":"IndexedGraph(A::AbstractMatrix)\n\nConstruct an IndexedGraph from symmetric adjacency matrix A.\n\n\n\n\n\n","category":"method"},{"location":"graph/","page":"IndexedGraph","title":"IndexedGraph","text":"inedges\noutedges\nget_edge(g::IndexedGraph, src::Integer, dst::Integer)","category":"page"},{"location":"graph/#IndexedGraphs.inedges","page":"IndexedGraph","title":"IndexedGraphs.inedges","text":"inedges(g::IndexedGraph, i::Integer)\n\nReturn a lazy iterators to the edges incident to i with i as the destination.\n\n\n\n\n\ninedges(g::AbstractIndexedBiDiGraph, i::Integer)\n\nReturn a lazy iterator to the edges ingoing to node i in g.\n\n\n\n\n\n","category":"function"},{"location":"graph/#IndexedGraphs.outedges","page":"IndexedGraph","title":"IndexedGraphs.outedges","text":"outedges(g::IndexedGraph, i::Integer)\n\nReturn a lazy iterators to the edges incident to i with i as the source.\n\n\n\n\n\noutedges(g::AbstractIndexedDiGraph, i::Integer)\n\nReturn a lazy iterator to the edges outgoing from node i in g.\n\n\n\n\n\n","category":"function"},{"location":"graph/#IndexedGraphs.get_edge-Tuple{IndexedGraph, Integer, Integer}","page":"IndexedGraph","title":"IndexedGraphs.get_edge","text":"get_edge(g::IndexedGraph, src::Integer, dst::Integer)\nget_edge(g::IndexedGraph, id::Integer)\n\nGet edge given source and destination or given edge index.\n\n\n\n\n\n","category":"method"},{"location":"graph/#Overrides-from-Graphs.jl","page":"IndexedGraph","title":"Overrides from Graphs.jl","text":"","category":"section"},{"location":"graph/","page":"IndexedGraph","title":"IndexedGraph","text":"edges(g::IndexedGraph, i::Integer)","category":"page"},{"location":"bidigraph/#IndexedBiDiGraph","page":"IndexedBiDiGraph","title":"IndexedBiDiGraph","text":"","category":"section"},{"location":"bidigraph/","page":"IndexedBiDiGraph","title":"IndexedBiDiGraph","text":"A type representing directed graphs. ","category":"page"},{"location":"bidigraph/","page":"IndexedBiDiGraph","title":"IndexedBiDiGraph","text":"Use this when you need to access both inedges and outedges (or inneighbors and outneighbors). For a lighter data structure check out IndexedDiGraph.","category":"page"},{"location":"bidigraph/","page":"IndexedBiDiGraph","title":"IndexedBiDiGraph","text":"IndexedBiDiGraph","category":"page"},{"location":"bidigraph/#IndexedGraphs.IndexedBiDiGraph","page":"IndexedBiDiGraph","title":"IndexedGraphs.IndexedBiDiGraph","text":"IndexedBiDiGraph{T<:Integer} <: AbstractIndexedDiGraph{T}\n\nA type representing a sparse directed graph with access to both outedges and inedges.\n\nFIELDS\n\nA – square matrix filled with NullNumbers. A[i,j] corresponds to edge j=>i.\nX – square matrix for efficient access by row. X[j,i] points to the index of element A[i,j] in A.nzval. \n\n\n\n\n\n","category":"type"},{"location":"bidigraph/","page":"IndexedBiDiGraph","title":"IndexedBiDiGraph","text":"IndexedBiDiGraph(A::AbstractMatrix)","category":"page"},{"location":"bidigraph/#IndexedGraphs.IndexedBiDiGraph-Tuple{AbstractMatrix}","page":"IndexedBiDiGraph","title":"IndexedGraphs.IndexedBiDiGraph","text":"IndexedBiDiGraph(A::AbstractMatrix)\n\nConstruct an IndexedBiDiGraph from the adjacency matrix A. \n\nIndexedBiDiGraph internally stores the transpose of A. To avoid overhead due to the transposition, use IndexedBiDiGraph(transpose(At)) where At is the  transpose of A.\n\n\n\n\n\n","category":"method"},{"location":"bidigraph/","page":"IndexedBiDiGraph","title":"IndexedBiDiGraph","text":"Example:","category":"page"},{"location":"bidigraph/","page":"IndexedBiDiGraph","title":"IndexedBiDiGraph","text":"using SparseArrays, IndexedGraphs\nAt = sprand(100, 100, 0.1)           # At[i,j] corresponds to edge j=>i\ng = IndexedBiDiGraph(transpose(At))  \ng.A.rowval === At.rowval","category":"page"},{"location":"digraph/#IndexedDiGraph","page":"IndexedDiGraph","title":"IndexedDiGraph","text":"","category":"section"},{"location":"digraph/","page":"IndexedDiGraph","title":"IndexedDiGraph","text":"A type representing directed graphs. ","category":"page"},{"location":"digraph/","page":"IndexedDiGraph","title":"IndexedDiGraph","text":"Use this when you need to access only outedges and outneighbors.  If you also need access to inedges and inneighbors, check out IndexedBiDiGraph.","category":"page"},{"location":"digraph/","page":"IndexedDiGraph","title":"IndexedDiGraph","text":"IndexedDiGraph","category":"page"},{"location":"digraph/#IndexedGraphs.IndexedDiGraph","page":"IndexedDiGraph","title":"IndexedGraphs.IndexedDiGraph","text":"IndexedDiGraph{T<:Integer} <: AbstractIndexedDiGraph{T}\n\nA type representing a sparse directed graph with access only to outedges.\n\nFIELDS\n\nA – square matrix filled with NullNumbers. A[i,j] corresponds to an edge j=>i\n\n\n\n\n\n","category":"type"},{"location":"digraph/","page":"IndexedDiGraph","title":"IndexedDiGraph","text":"IndexedDiGraph(A::AbstractMatrix)","category":"page"},{"location":"digraph/#IndexedGraphs.IndexedDiGraph-Tuple{AbstractMatrix}","page":"IndexedDiGraph","title":"IndexedGraphs.IndexedDiGraph","text":"IndexedDiGraph(A::AbstractMatrix)\n\nConstructs a IndexedDiGraph from the adjacency matrix A.\n\nIndexedDiGraph internally stores the transpose of A. To avoid overhead due to the transposition, use IndexedDiGraph(transpose(At)) where At is the  transpose of A.\n\n\n\n\n\n","category":"method"},{"location":"factor/#FactorGraph","page":"FactorGraph","title":"FactorGraph","text":"","category":"section"},{"location":"factor/","page":"FactorGraph","title":"FactorGraph","text":"A type representing Factor Graphs.","category":"page"},{"location":"factor/","page":"FactorGraph","title":"FactorGraph","text":"FactorGraph","category":"page"},{"location":"factor/#IndexedGraphs.FactorGraph","page":"FactorGraph","title":"IndexedGraphs.FactorGraph","text":"FactorGraph{T<:Integer} <: AbstractIndexedGraph{T}\n\nA type representing a sparse factor graph.\n\nFIELDS\n\nA – adjacency matrix filled with NullNumbers. Rows are factors, columns are variables\nX – square matrix for efficient access by row. X[j,i] points to the index of element A[i,j] in A.nzval. \n\n\n\n\n\n","category":"type"},{"location":"factor/","page":"FactorGraph","title":"FactorGraph","text":"FactorGraph(A::AbstractMatrix)","category":"page"},{"location":"factor/#IndexedGraphs.FactorGraph-Tuple{AbstractMatrix}","page":"FactorGraph","title":"IndexedGraphs.FactorGraph","text":"FactorGraph(A::AbstractMatrix)\n\nConstruct a factor graph from adjacency matrix A with the convention that rows are factors, columns are variables\n\n\n\n\n\n","category":"method"},{"location":"factor/","page":"FactorGraph","title":"FactorGraph","text":"using IndexedGraphs.VariableOrFactor","category":"page"},{"location":"factor/","page":"FactorGraph","title":"FactorGraph","text":"get_edge(g::FactorGraph, x::VariableOrFactor, y::VariableOrFactor)\nbipartite_view(g::FactorGraph, T::DataType=Int)","category":"page"},{"location":"factor/#IndexedGraphs.get_edge-Tuple{FactorGraph, VariableOrFactor, VariableOrFactor}","page":"FactorGraph","title":"IndexedGraphs.get_edge","text":"get_edge(g::FactorGraph, v::Variable, f::Factor)\nget_edge(g::FactorGraph, f::Factor, v::Variable)\nget_edge(g::FactorGraph, id::Integer)\n\nGet edge given variable and factor or given edge index\n\n\n\n\n\n","category":"method"},{"location":"factor/#IndexedGraphs.bipartite_view","page":"FactorGraph","title":"IndexedGraphs.bipartite_view","text":"bipartite_view(g::FactorGraph, T::DataType=Int)\n\nConstruct an undirected bipartite IndexedGraph where vertices are the concatenation of the variables and factors of the original FactorGraph.  Edge indices are preserved in the bottom-left block.\n\n\n\n\n\n","category":"function"},{"location":"factor/#Overrides-from-Graphs.jl","page":"FactorGraph","title":"Overrides from Graphs.jl","text":"","category":"section"},{"location":"factor/","page":"FactorGraph","title":"FactorGraph","text":"neighbors(g::FactorGraph, v::Variable)\nneighbors(g::FactorGraph, f::Factor)\nedges(g::FactorGraph, v::Variable)\nedges(g::FactorGraph, f::Factor)\nadjacency_matrix(g::FactorGraph, T::DataType=Int)","category":"page"},{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [IndexedGraphs]","category":"page"},{"location":"reference/#IndexedGraphs.AbstractIndexedDiGraph","page":"Reference","title":"IndexedGraphs.AbstractIndexedDiGraph","text":"AbstractIndexedDiGraph{T} <: AbstractIndexedGraph{T}\n\nAbstract type for representing directed graphs.\n\n\n\n\n\n","category":"type"},{"location":"reference/#IndexedGraphs.AbstractIndexedEdge","page":"Reference","title":"IndexedGraphs.AbstractIndexedEdge","text":"AbstractIndexedEdge{T<:Integer} <: AbstractEdge{T}\n\nAbstract type for indexed edge. AbstractIndexedEdge{T}s must have the following elements:\n\nidx::T integer positive index \n\n\n\n\n\n","category":"type"},{"location":"reference/#IndexedGraphs.AbstractIndexedGraph","page":"Reference","title":"IndexedGraphs.AbstractIndexedGraph","text":"AbstractIndexedGraph{T} <: AbstractGraph{T}\n\nAn abstract type representing an indexed graph. AbstractIndexedGraphs must have the following elements:\n\nA::SparseMatrixCSC adjacency matrix\n\n\n\n\n\n","category":"type"},{"location":"reference/#IndexedGraphs.Factor","page":"Reference","title":"IndexedGraphs.Factor","text":"Factor{T<:Integer}\n\nWraps an index to specify that it is the index of a factor.  See e.g. neighbors(g::FactorGraph, f::Factor)\n\n\n\n\n\n","category":"type"},{"location":"reference/#IndexedGraphs.FactorGraphEdge","page":"Reference","title":"IndexedGraphs.FactorGraphEdge","text":"FactorGraphEdge{T<:Integer} <: AbstractIndexedEdge{T}\n\nEdge type for FactorGraphs.\n\nSafe constructors\n\nFactorGraphEdge(f::Factor, v::Variable, idx::Integer)\nFactorGraphEdge(v::Variable, f::Factor, idx::Integer)\n\nExample:\n\n    FactorGraphEdge(Factor(2), Variable(4), 3)\n\n\n\n\n\n","category":"type"},{"location":"reference/#IndexedGraphs.IndexedEdge","page":"Reference","title":"IndexedGraphs.IndexedEdge","text":"IndexedEdge{T<:Integer} <: AbstractIndexedEdge{T}\n\nEdge type for IndexedGraphs. Edge indices can be used to access edge  properties stored in separate containers.\n\n\n\n\n\n","category":"type"},{"location":"reference/#IndexedGraphs.Variable","page":"Reference","title":"IndexedGraphs.Variable","text":"Variable{T<:Integer}\n\nWraps an index to specify that it is the index of a variable.  See e.g. neighbors(g::FactorGraph, v::Variable)\n\n\n\n\n\n","category":"type"},{"location":"reference/#Graphs.LinAlg.adjacency_matrix","page":"Reference","title":"Graphs.LinAlg.adjacency_matrix","text":"adjacency_matrix(g::FactorGraph, T::DataType=Int)\n\nReturn the symmetric adjacency matrix of size nvariables(g) + nfactors(g)  where no distinction is made between variable and factor nodes. Edge indices are preserved in the bottom-left block.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphs.edges-Tuple{FactorGraph, Factor}","page":"Reference","title":"Graphs.edges","text":"neighbors(g::FactorGraph, f::Factor)\n\nReturn a lazy iterator to the edges incident to factor f\n\n\n\n\n\n","category":"method"},{"location":"reference/#Graphs.edges-Tuple{FactorGraph, Variable}","page":"Reference","title":"Graphs.edges","text":"edges(g::FactorGraph, v::Variable)\n\nReturn a lazy iterator to the edges incident to variable v\n\n\n\n\n\n","category":"method"},{"location":"reference/#Graphs.edges-Tuple{IndexedGraph, Integer}","page":"Reference","title":"Graphs.edges","text":"edges(g::IndexedGraph, i::Integer)\n\nReturn a lazy iterators to the edges incident to i.\n\nBy default unordered edges sort source and destination nodes in increasing order. See outedges and inedges if you need otherwise.\n\n\n\n\n\n","category":"method"},{"location":"reference/#Graphs.neighbors-Tuple{FactorGraph, Factor}","page":"Reference","title":"Graphs.neighbors","text":"neighbors(g::FactorGraph, f::Factor)\n\nReturn a lazy iterator to the neighbors of factor f\n\n\n\n\n\n","category":"method"},{"location":"reference/#Graphs.neighbors-Tuple{FactorGraph, Variable}","page":"Reference","title":"Graphs.neighbors","text":"neighbors(g::FactorGraph, v::Variable)\n\nReturn a lazy iterator to the neighbors of variable v\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.get_edge-Tuple{AbstractIndexedDiGraph, Integer, Integer}","page":"Reference","title":"IndexedGraphs.get_edge","text":"get_edge(g::AbstractIndexedDiGraph, src::Integer, dst::Integer)\nget_edge(g::AbstractIndexedDiGraph, id::Integer)\n\nGet edge given source and destination or given edge index.\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.inedges-Tuple{IndexedBiDiGraph, Integer}","page":"Reference","title":"IndexedGraphs.inedges","text":"inedges(g::AbstractIndexedBiDiGraph, i::Integer)\n\nReturn a lazy iterator to the edges ingoing to node i in g.\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.inedges-Tuple{IndexedGraph, Integer}","page":"Reference","title":"IndexedGraphs.inedges","text":"inedges(g::IndexedGraph, i::Integer)\n\nReturn a lazy iterators to the edges incident to i with i as the destination.\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.outedges-Tuple{AbstractIndexedDiGraph, Integer}","page":"Reference","title":"IndexedGraphs.outedges","text":"outedges(g::AbstractIndexedDiGraph, i::Integer)\n\nReturn a lazy iterator to the edges outgoing from node i in g.\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.outedges-Tuple{IndexedGraph, Integer}","page":"Reference","title":"IndexedGraphs.outedges","text":"outedges(g::IndexedGraph, i::Integer)\n\nReturn a lazy iterators to the edges incident to i with i as the source.\n\n\n\n\n\n","category":"method"},{"location":"#IndexedGraphs.jl","page":"Home","title":"IndexedGraphs.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package defines three basic types of Graphs:","category":"page"},{"location":"","page":"Home","title":"Home","text":"IndexedGraph\nIndexedDiGraph\nIndexedBidiGraph","category":"page"},{"location":"","page":"Home","title":"Home","text":"In addition, it provides a FactorGraph type.","category":"page"},{"location":"","page":"Home","title":"Home","text":"They all comply with the Developing Alternate Graph Types rules for subtyping from Graphs.AbstractGraph.","category":"page"}]
}
