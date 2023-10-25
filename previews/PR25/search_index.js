var documenterSearchIndex = {"docs":
[{"location":"graph/#IndexedGraph","page":"IndexedGraph","title":"IndexedGraph","text":"","category":"section"},{"location":"graph/","page":"IndexedGraph","title":"IndexedGraph","text":"A type representing undirected graphs.","category":"page"},{"location":"graph/","page":"IndexedGraph","title":"IndexedGraph","text":"IndexedGraph","category":"page"},{"location":"graph/#IndexedGraphs.IndexedGraph","page":"IndexedGraph","title":"IndexedGraphs.IndexedGraph","text":"IndexedGraph{T<:Integer} <: AbstractIndexedGraph{T}\n\nA type representing a sparse undirected graph.\n\nFIELDS\n\nA – square adjacency matrix. A[i,j] == A[j,i] contains the unique index associated to unidrected edge (i,j)\n\n\n\n\n\n","category":"type"},{"location":"graph/","page":"IndexedGraph","title":"IndexedGraph","text":"IndexedGraph(A::AbstractMatrix)","category":"page"},{"location":"graph/#IndexedGraphs.IndexedGraph-Tuple{AbstractMatrix}","page":"IndexedGraph","title":"IndexedGraphs.IndexedGraph","text":"IndexedGraph(A::AbstractMatrix)\n\nConstruct an IndexedGraph from symmetric adjacency matrix A.\n\n\n\n\n\n","category":"method"},{"location":"graph/","page":"IndexedGraph","title":"IndexedGraph","text":"inedges(g::IndexedGraph, i::Integer)\noutedges(g::IndexedGraph, i::Integer)\nget_edge(g::IndexedGraph, src::Integer, dst::Integer)","category":"page"},{"location":"graph/#IndexedGraphs.inedges-Tuple{IndexedGraph, Integer}","page":"IndexedGraph","title":"IndexedGraphs.inedges","text":"inedges(g::IndexedGraph, i::Integer)\n\nReturn a lazy iterators to the edges incident to i with i as the destination.\n\n\n\n\n\n","category":"method"},{"location":"graph/#IndexedGraphs.outedges-Tuple{IndexedGraph, Integer}","page":"IndexedGraph","title":"IndexedGraphs.outedges","text":"outedges(g::IndexedGraph, i::Integer)\n\nReturn a lazy iterators to the edges incident to i with i as the source.\n\n\n\n\n\n","category":"method"},{"location":"graph/#IndexedGraphs.get_edge-Tuple{IndexedGraph, Integer, Integer}","page":"IndexedGraph","title":"IndexedGraphs.get_edge","text":"get_edge(g::IndexedGraph, src::Integer, dst::Integer)\nget_edge(g::IndexedGraph, id::Integer)\n\nGet edge given source and destination or given edge index.\n\n\n\n\n\n","category":"method"},{"location":"graph/#Overrides-from-Graphs.jl","page":"IndexedGraph","title":"Overrides from Graphs.jl","text":"","category":"section"},{"location":"graph/","page":"IndexedGraph","title":"IndexedGraph","text":"edges(g::IndexedGraph, i::Integer)","category":"page"},{"location":"graph/#Graphs.edges-Tuple{IndexedGraph, Integer}","page":"IndexedGraph","title":"Graphs.edges","text":"edges(g::IndexedGraph, i::Integer)\n\nReturn a lazy iterators to the edges incident to i.\n\nBy default unordered edges sort source and destination nodes in increasing order. See outedges and inedges if you need otherwise.\n\n\n\n\n\n","category":"method"},{"location":"bidigraph/#IndexedBiDiGraph","page":"IndexedBiDiGraph","title":"IndexedBiDiGraph","text":"","category":"section"},{"location":"bidigraph/","page":"IndexedBiDiGraph","title":"IndexedBiDiGraph","text":"A type representing directed graphs. ","category":"page"},{"location":"bidigraph/","page":"IndexedBiDiGraph","title":"IndexedBiDiGraph","text":"Use this when you need to access both inedges and outedges (or inneighbors and outneighbors). For a lighter data structure check out IndexedDiGraph.","category":"page"},{"location":"bidigraph/","page":"IndexedBiDiGraph","title":"IndexedBiDiGraph","text":"IndexedBiDiGraph","category":"page"},{"location":"bidigraph/#IndexedGraphs.IndexedBiDiGraph","page":"IndexedBiDiGraph","title":"IndexedGraphs.IndexedBiDiGraph","text":"IndexedBiDiGraph{T<:Integer} <: AbstractIndexedDiGraph{T}\n\nA type representing a sparse directed graph with access to both outedges and inedges.\n\nFIELDS\n\nA – square matrix filled with NullNumbers. A[i,j] corresponds to edge j=>i.\nX – square matrix for efficient access by row. X[j,i] points to the index of element A[i,j] in A.nzval. \n\n\n\n\n\n","category":"type"},{"location":"bidigraph/","page":"IndexedBiDiGraph","title":"IndexedBiDiGraph","text":"IndexedBiDiGraph(A::AbstractMatrix)","category":"page"},{"location":"bidigraph/#IndexedGraphs.IndexedBiDiGraph-Tuple{AbstractMatrix}","page":"IndexedBiDiGraph","title":"IndexedGraphs.IndexedBiDiGraph","text":"IndexedBiDiGraph(A::AbstractMatrix)\n\nConstruct an IndexedBiDiGraph from the adjacency matrix A. \n\nIndexedBiDiGraph internally stores the transpose of A. To avoid overhead due to the transposition, use IndexedBiDiGraph(transpose(At)) where At is the  transpose of A.\n\n\n\n\n\n","category":"method"},{"location":"bidigraph/","page":"IndexedBiDiGraph","title":"IndexedBiDiGraph","text":"Example:","category":"page"},{"location":"bidigraph/","page":"IndexedBiDiGraph","title":"IndexedBiDiGraph","text":"using SparseArrays, IndexedGraphs\nAt = sprand(100, 100, 0.1)           # At[i,j] corresponds to edge j=>i\nfor i in 1:100; At[i,i] = 0; end\ndropzeros!(At)\ng = IndexedBiDiGraph(transpose(At))  \ng.A.rowval === At.rowval","category":"page"},{"location":"digraph/#IndexedDiGraph","page":"IndexedDiGraph","title":"IndexedDiGraph","text":"","category":"section"},{"location":"digraph/","page":"IndexedDiGraph","title":"IndexedDiGraph","text":"A type representing directed graphs. ","category":"page"},{"location":"digraph/","page":"IndexedDiGraph","title":"IndexedDiGraph","text":"Use this when you need to access only outedges and outneighbors.  If you also need access to inedges and inneighbors, check out IndexedBiDiGraph.","category":"page"},{"location":"digraph/","page":"IndexedDiGraph","title":"IndexedDiGraph","text":"IndexedDiGraph","category":"page"},{"location":"digraph/#IndexedGraphs.IndexedDiGraph","page":"IndexedDiGraph","title":"IndexedGraphs.IndexedDiGraph","text":"IndexedDiGraph{T<:Integer} <: AbstractIndexedDiGraph{T}\n\nA type representing a sparse directed graph with access only to outedges.\n\nFIELDS\n\nA – square matrix filled with NullNumbers. A[i,j] corresponds to an edge j=>i\n\n\n\n\n\n","category":"type"},{"location":"digraph/","page":"IndexedDiGraph","title":"IndexedDiGraph","text":"IndexedDiGraph(A::AbstractMatrix)","category":"page"},{"location":"digraph/#IndexedGraphs.IndexedDiGraph-Tuple{AbstractMatrix}","page":"IndexedDiGraph","title":"IndexedGraphs.IndexedDiGraph","text":"IndexedDiGraph(A::AbstractMatrix)\n\nConstructs a IndexedDiGraph from the adjacency matrix A.\n\nIndexedDiGraph internally stores the transpose of A. To avoid overhead due to the transposition, use IndexedDiGraph(transpose(At)) where At is the  transpose of A.\n\n\n\n\n\n","category":"method"},{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [IndexedGraphs]","category":"page"},{"location":"reference/#IndexedGraphs.AbstractIndexedDiGraph","page":"Reference","title":"IndexedGraphs.AbstractIndexedDiGraph","text":"AbstractIndexedDiGraph{T} <: AbstractIndexedGraph{T}\n\nAbstract type for representing directed graphs.\n\n\n\n\n\n","category":"type"},{"location":"reference/#IndexedGraphs.AbstractIndexedEdge","page":"Reference","title":"IndexedGraphs.AbstractIndexedEdge","text":"AbstractIndexedEdge{T<:Integer} <: AbstractEdge{T}\n\nAbstract type for indexed edge. AbstractIndexedEdge{T}s must have the following elements:\n\nidx::T integer positive index \n\n\n\n\n\n","category":"type"},{"location":"reference/#IndexedGraphs.AbstractIndexedGraph","page":"Reference","title":"IndexedGraphs.AbstractIndexedGraph","text":"AbstractIndexedGraph{T} <: AbstractGraph{T}\n\nAn abstract type representing an indexed graph. AbstractIndexedGraphs must have the following elements:\n\nA::SparseMatrixCSC adjacency matrix\n\n\n\n\n\n","category":"type"},{"location":"reference/#IndexedGraphs.IndexedBiDiGraph-Tuple{AbstractSimpleGraph}","page":"Reference","title":"IndexedGraphs.IndexedBiDiGraph","text":"IndexedBiDiGraph(A::AbstractSimpleGraph)\n\nConstruct an IndexedBiDiGraph from any AbstractSimpleGraph (Graphs.jl),  directed or otherwise.\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.IndexedDiGraph-Tuple{AbstractSimpleGraph}","page":"Reference","title":"IndexedGraphs.IndexedDiGraph","text":"IndexedDiGraph(A::AbstractSimpleGraph)\n\nConstruct an IndexedDiGraph from any AbstractSimpleGraph (Graphs.jl),  directed or otherwise.\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.IndexedEdge","page":"Reference","title":"IndexedGraphs.IndexedEdge","text":"IndexedEdge{T<:Integer} <: AbstractIndexedEdge{T}\n\nEdge type for IndexedGraphs. Edge indices can be used to access edge  properties stored in separate containers.\n\n\n\n\n\n","category":"type"},{"location":"reference/#IndexedGraphs.IndexedGraph-Tuple{SimpleGraph}","page":"Reference","title":"IndexedGraphs.IndexedGraph","text":"IndexedGraph(A::SimpleGraph)\n\nConstruct an IndexedGraph from undirected SimpleGraph (Graphs.jl).\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.get_edge-Tuple{AbstractIndexedDiGraph, Integer, Integer}","page":"Reference","title":"IndexedGraphs.get_edge","text":"get_edge(g::AbstractIndexedDiGraph, src::Integer, dst::Integer)\nget_edge(g::AbstractIndexedDiGraph, id::Integer)\n\nGet edge given source and destination or given edge index.\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.inedges-Tuple{IndexedBiDiGraph, Integer}","page":"Reference","title":"IndexedGraphs.inedges","text":"inedges(g::AbstractIndexedBiDiGraph, i::Integer)\n\nReturn a lazy iterator to the edges ingoing to node i in g.\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.linearindex-Tuple{BipartiteIndexedGraph, BipartiteGraphVertex{Right}}","page":"Reference","title":"IndexedGraphs.linearindex","text":"linearindex(g::BipartiteIndexedGraph, v::BipartiteGraphVertex{LR}) where {LR<:LeftorRight}\nlinearindex(g::BipartiteIndexedGraph, i::Integer, ::Type{LR}) where LR<:LeftorRight\n\nReturn the linear index of a vertex, specified either by a BipartiteGraphVertex or by its index and block.\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.nv_left-Tuple{BipartiteIndexedGraph}","page":"Reference","title":"IndexedGraphs.nv_left","text":"nv_left(g::BipartiteIndexedGraph)\n\nReturn the number of vertices in the left block\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.nv_right-Tuple{BipartiteIndexedGraph}","page":"Reference","title":"IndexedGraphs.nv_right","text":"nv_right(g::BipartiteIndexedGraph)\n\nReturn the number of vertices in the right block\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.outedges-Tuple{AbstractIndexedDiGraph, Integer}","page":"Reference","title":"IndexedGraphs.outedges","text":"outedges(g::AbstractIndexedDiGraph, i::Integer)\n\nReturn a lazy iterator to the edges outgoing from node i in g.\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.vertex-Union{Tuple{LR}, Tuple{Integer, Type{LR}}} where LR<:Union{Left, Right}","page":"Reference","title":"IndexedGraphs.vertex","text":"vertex(i::Integer, ::Type{<:LeftorRight})\n\nBuild a BipartiteGraphVertex\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.vertex_left-Tuple{BipartiteIndexedGraph, IndexedGraphs.IndexedEdge}","page":"Reference","title":"IndexedGraphs.vertex_left","text":"vertex_left(g::BipartiteIndexedGraph, e::IndexedEdge)\n\nReturn the (in-block) index of the left-block vertex in e.\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.vertex_right-Tuple{BipartiteIndexedGraph, IndexedGraphs.IndexedEdge}","page":"Reference","title":"IndexedGraphs.vertex_right","text":"vertex_right(g::BipartiteIndexedGraph, e::IndexedEdge)\n\nReturn the (in-block) index of the right-block vertex in e.\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.vertices_left-Tuple{BipartiteIndexedGraph}","page":"Reference","title":"IndexedGraphs.vertices_left","text":"vertices_left(g::BipartiteIndexedGraph)\n\nReturn a lazy iterator to the vertices in the left block\n\n\n\n\n\n","category":"method"},{"location":"reference/#IndexedGraphs.vertices_right-Tuple{BipartiteIndexedGraph}","page":"Reference","title":"IndexedGraphs.vertices_right","text":"vertices right(g::BipartiteIndexedGraph)\n\nReturn a lazy iterator to the vertices in the right block\n\n\n\n\n\n","category":"method"},{"location":"bipartite/#BipartiteIndexedGraph","page":"BipartiteIndexedGraph","title":"BipartiteIndexedGraph","text":"","category":"section"},{"location":"bipartite/","page":"BipartiteIndexedGraph","title":"BipartiteIndexedGraph","text":"A graph is bipartite if the set of vertices can be partitioned into two blocks such that there is no edge between vertices of different blocks. Here we adopt the notation of referring to these as the \"left\" and \"right\" block.","category":"page"},{"location":"bipartite/","page":"BipartiteIndexedGraph","title":"BipartiteIndexedGraph","text":"A BipartiteIndexedGraph behaves just like a bipartite, undirected IndexedGraph, with two differences:","category":"page"},{"location":"bipartite/","page":"BipartiteIndexedGraph","title":"BipartiteIndexedGraph","text":"Vertices can be indexed as usual via their integer index (which is called here a linearindex), or via a BipartiteGraphVertex, i.e. by specifying Left or Right and the integer index of that vertex within its block. The typical use case is that where one has two vectors storing vertex properties of the two blocks, possibly with different eltypes, and each with indices starting at one.\nThe adjacency matrix of a bipartite graph (possibly after permutation of the vertex indices) is made of two specular rectangular submatrices A. Only one of these is stored, leading to a slight improvement in efficiency.","category":"page"},{"location":"bipartite/","page":"BipartiteIndexedGraph","title":"BipartiteIndexedGraph","text":"BipartiteIndexedGraphs use the same edge type IndexedEdge as the other AbstractIndexedGraphs, which stores source and destination as linear indices. To retrieve the in-block indices of the two vertices incident on an edge, use vertex_left, vertex_right.","category":"page"},{"location":"bipartite/","page":"BipartiteIndexedGraph","title":"BipartiteIndexedGraph","text":"BipartiteIndexedGraph","category":"page"},{"location":"bipartite/#IndexedGraphs.BipartiteIndexedGraph","page":"BipartiteIndexedGraph","title":"IndexedGraphs.BipartiteIndexedGraph","text":"BipartiteIndexedGraph{T<:Integer} <: AbstractIndexedGraph{T}\n\nA type representing a sparse, undirected bipartite graph.\n\nFIELDS\n\nA – adjacency matrix filled with NullNumbers. Rows are vertices belonging to the left block, columns to the right block\nX – square matrix for efficient access by row. X[j,i] points to the index of element A[i,j] in A.nzval. \n\n\n\n\n\n","category":"type"},{"location":"bipartite/","page":"BipartiteIndexedGraph","title":"BipartiteIndexedGraph","text":"BipartiteIndexedGraph(A::AbstractMatrix)\nBipartiteIndexedGraph(g::AbstractGraph)\nnv_left\nnv_right\nLeft\nRight\nLeftorRight\nBipartiteGraphVertex\nvertex(i::Integer, ::Type{<:LeftorRight})\nvertex(g::BipartiteIndexedGraph, i::Integer)\nlinearindex\nvertices_left\nvertices_right\nvertex_left\nvertex_right\ninedges(g::BipartiteIndexedGraph, i::Integer)\ninedges(g::BipartiteIndexedGraph, l::BipartiteGraphVertex{Left})\noutedges(g::BipartiteIndexedGraph, i::Integer)\noutedges(g::BipartiteIndexedGraph, l::BipartiteGraphVertex{Left})","category":"page"},{"location":"bipartite/#IndexedGraphs.BipartiteIndexedGraph-Tuple{AbstractMatrix}","page":"BipartiteIndexedGraph","title":"IndexedGraphs.BipartiteIndexedGraph","text":"BipartiteIndexedGraph(A::AbstractMatrix)\n\nConstruct a BipartiteIndexedGraph from adjacency matrix A with the convention that rows are vertices belonging to the left block, columns to the right block\n\n\n\n\n\n","category":"method"},{"location":"bipartite/#IndexedGraphs.BipartiteIndexedGraph-Tuple{AbstractGraph}","page":"BipartiteIndexedGraph","title":"IndexedGraphs.BipartiteIndexedGraph","text":"BipartiteIndexedGraph(g::AbstractGraph)\n\nBuild a BipartiteIndexedGraph from any undirected, bipartite graph.\n\n\n\n\n\n","category":"method"},{"location":"bipartite/#IndexedGraphs.nv_left","page":"BipartiteIndexedGraph","title":"IndexedGraphs.nv_left","text":"nv_left(g::BipartiteIndexedGraph)\n\nReturn the number of vertices in the left block\n\n\n\n\n\n","category":"function"},{"location":"bipartite/#IndexedGraphs.nv_right","page":"BipartiteIndexedGraph","title":"IndexedGraphs.nv_right","text":"nv_right(g::BipartiteIndexedGraph)\n\nReturn the number of vertices in the right block\n\n\n\n\n\n","category":"function"},{"location":"bipartite/#IndexedGraphs.Left","page":"BipartiteIndexedGraph","title":"IndexedGraphs.Left","text":"Left\n\nSingleton type used to represent a vertex belonging to the left block in a BipartiteGraphVertex\n\n\n\n\n\n","category":"type"},{"location":"bipartite/#IndexedGraphs.Right","page":"BipartiteIndexedGraph","title":"IndexedGraphs.Right","text":"Right\n\nSingleton type used to represent a vertex belonging to the Right block in a BipartiteGraphVertex\n\n\n\n\n\n","category":"type"},{"location":"bipartite/#IndexedGraphs.LeftorRight","page":"BipartiteIndexedGraph","title":"IndexedGraphs.LeftorRight","text":"LeftorRight\n\nLeftorRight = Union{Left, Right}\n\n\n\n\n\n","category":"type"},{"location":"bipartite/#IndexedGraphs.BipartiteGraphVertex","page":"BipartiteIndexedGraph","title":"IndexedGraphs.BipartiteGraphVertex","text":"BipartiteGraphVertex\n\nA BipartiteGraphVertex{LR<:LeftorRight,T<:Integer} represents a vertex in a bipartite graph.\n\nPARAMETERS\n\nLR – Either Left or Right\n\nFIELDS\n\ni – The index of the vertex within its block.\n\n\n\n\n\n","category":"type"},{"location":"bipartite/#IndexedGraphs.vertex-Tuple{Integer, Type{<:Union{Left, Right}}}","page":"BipartiteIndexedGraph","title":"IndexedGraphs.vertex","text":"vertex(i::Integer, ::Type{<:LeftorRight})\n\nBuild a BipartiteGraphVertex\n\n\n\n\n\n","category":"method"},{"location":"bipartite/#IndexedGraphs.vertex-Tuple{BipartiteIndexedGraph, Integer}","page":"BipartiteIndexedGraph","title":"IndexedGraphs.vertex","text":"vertex(g::BipartiteIndexedGraph, i::Integer)\n\nBuild the BipartiteGraphVertex corresponding to linear index i. Throws an error if i is not in the range of vertices of g\n\n\n\n\n\n","category":"method"},{"location":"bipartite/#IndexedGraphs.linearindex","page":"BipartiteIndexedGraph","title":"IndexedGraphs.linearindex","text":"linearindex(g::BipartiteIndexedGraph, v::BipartiteGraphVertex{LR}) where {LR<:LeftorRight}\nlinearindex(g::BipartiteIndexedGraph, i::Integer, ::Type{LR}) where LR<:LeftorRight\n\nReturn the linear index of a vertex, specified either by a BipartiteGraphVertex or by its index and block.\n\n\n\n\n\n","category":"function"},{"location":"bipartite/#IndexedGraphs.vertices_left","page":"BipartiteIndexedGraph","title":"IndexedGraphs.vertices_left","text":"vertices_left(g::BipartiteIndexedGraph)\n\nReturn a lazy iterator to the vertices in the left block\n\n\n\n\n\n","category":"function"},{"location":"bipartite/#IndexedGraphs.vertices_right","page":"BipartiteIndexedGraph","title":"IndexedGraphs.vertices_right","text":"vertices right(g::BipartiteIndexedGraph)\n\nReturn a lazy iterator to the vertices in the right block\n\n\n\n\n\n","category":"function"},{"location":"bipartite/#IndexedGraphs.vertex_left","page":"BipartiteIndexedGraph","title":"IndexedGraphs.vertex_left","text":"vertex_left(g::BipartiteIndexedGraph, e::IndexedEdge)\n\nReturn the (in-block) index of the left-block vertex in e.\n\n\n\n\n\n","category":"function"},{"location":"bipartite/#IndexedGraphs.vertex_right","page":"BipartiteIndexedGraph","title":"IndexedGraphs.vertex_right","text":"vertex_right(g::BipartiteIndexedGraph, e::IndexedEdge)\n\nReturn the (in-block) index of the right-block vertex in e.\n\n\n\n\n\n","category":"function"},{"location":"bipartite/#IndexedGraphs.inedges-Tuple{BipartiteIndexedGraph, Integer}","page":"BipartiteIndexedGraph","title":"IndexedGraphs.inedges","text":"inedges(g::BipartiteIndexedGraph, i::Integer)\n\nReturn a lazy iterator to the edges incident on a variable i specified by its linear index, with i as the destination.\n\n\n\n\n\n","category":"method"},{"location":"bipartite/#IndexedGraphs.inedges-Tuple{BipartiteIndexedGraph, BipartiteGraphVertex{Left}}","page":"BipartiteIndexedGraph","title":"IndexedGraphs.inedges","text":"inedges(g::BipartiteIndexedGraph, v::BipartiteGraphVertex{<:LeftorRight})\n\nReturn a lazy iterator to the edges incident on a variable v specified by a BipartiteGraphVertex, with v as the destination. \n\n\n\n\n\n","category":"method"},{"location":"bipartite/#IndexedGraphs.outedges-Tuple{BipartiteIndexedGraph, Integer}","page":"BipartiteIndexedGraph","title":"IndexedGraphs.outedges","text":"outedges(g::BipartiteIndexedGraph, i::Integer)\n\nReturn a lazy iterator to the edges incident on a variable i specified by its linear index, with i as the source.\n\n\n\n\n\n","category":"method"},{"location":"bipartite/#IndexedGraphs.outedges-Tuple{BipartiteIndexedGraph, BipartiteGraphVertex{Left}}","page":"BipartiteIndexedGraph","title":"IndexedGraphs.outedges","text":"outedges(g::BipartiteIndexedGraph, v::BipartiteGraphVertex{<:LeftorRight})\n\nReturn a lazy iterator to the edges incident on a variable v specified by a BipartiteGraphVertex, with v as the source. \n\n\n\n\n\n","category":"method"},{"location":"bipartite/#Overrides-from-Graphs.jl","page":"BipartiteIndexedGraph","title":"Overrides from Graphs.jl","text":"","category":"section"},{"location":"bipartite/","page":"BipartiteIndexedGraph","title":"BipartiteIndexedGraph","text":"inneighbors(g::BipartiteIndexedGraph, i::Integer)\ninneighbors(g::BipartiteIndexedGraph, l::BipartiteGraphVertex{Left})\noutneighbors(g::BipartiteIndexedGraph, i::Integer)\noutneighbors(g::BipartiteIndexedGraph, v::BipartiteGraphVertex)\nadjacency_matrix(g::BipartiteIndexedGraph, T::DataType=Int)","category":"page"},{"location":"bipartite/#Graphs.inneighbors-Tuple{BipartiteIndexedGraph, Integer}","page":"BipartiteIndexedGraph","title":"Graphs.inneighbors","text":"inneighbors(g::BipartiteIndexedGraph, i::Integer)\n\nReturn a lazy iterator to the neighbors of variable i specified by its linear index. \n\n\n\n\n\n","category":"method"},{"location":"bipartite/#Graphs.inneighbors-Tuple{BipartiteIndexedGraph, BipartiteGraphVertex{Left}}","page":"BipartiteIndexedGraph","title":"Graphs.inneighbors","text":"inneighbors(g::BipartiteIndexedGraph, v::BipartiteGraphVertex{<:LeftorRight})\n\nReturn a lazy iterator to the neighbors of variable v specified by a BipartiteGraphVertex. \n\n\n\n\n\n","category":"method"},{"location":"bipartite/#Graphs.outneighbors-Tuple{BipartiteIndexedGraph, Integer}","page":"BipartiteIndexedGraph","title":"Graphs.outneighbors","text":"outneighbors(g::BipartiteIndexedGraph, i::Integer)\n\nReturn a lazy iterator to the neighbors of variable i specified by its linear index. \n\n\n\n\n\n","category":"method"},{"location":"bipartite/#Graphs.outneighbors-Tuple{BipartiteIndexedGraph, BipartiteGraphVertex}","page":"BipartiteIndexedGraph","title":"Graphs.outneighbors","text":"outneighbors(g::BipartiteIndexedGraph, v::BipartiteGraphVertex{<:LeftorRight})\n\nReturn a lazy iterator to the neighbors of variable v specified by a BipartiteGraphVertex. \n\n\n\n\n\n","category":"method"},{"location":"bipartite/#Graphs.LinAlg.adjacency_matrix","page":"BipartiteIndexedGraph","title":"Graphs.LinAlg.adjacency_matrix","text":"adjacency_matrix(g::BipartiteIndexedGraph, T::DataType=Int)\n\nReturn the symmetric adjacency matrix of size nv(g) = nv_left(g) + nv_right(g)  where no distinction is made between left and right nodes.\n\n\n\n\n\n","category":"function"},{"location":"#IndexedGraphs.jl","page":"Home","title":"IndexedGraphs.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package defines three basic types of Graphs:","category":"page"},{"location":"","page":"Home","title":"Home","text":"IndexedGraph\nIndexedDiGraph\nIndexedBiDiGraph","category":"page"},{"location":"","page":"Home","title":"Home","text":"In addition, it provides a BipartiteIndexedGraph type.","category":"page"},{"location":"","page":"Home","title":"Home","text":"They all comply with the Developing Alternate Graph Types rules for subtyping from Graphs.AbstractGraph.","category":"page"}]
}
