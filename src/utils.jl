function _checksquare(A::AbstractMatrix)
    size(A,1) != size(A,2) && throw(ArgumentError("Matrix should be square"))
end

function _check_selfloops(A::AbstractMatrix)
    s = min(size(A)...)
    any(!iszero, A[i,i] for i=1:s) && throw(ArgumentError("Self loops are not allowed"))
end

# vector that maps edge i->j to edge j->i in A.nzval
function inverse_edge_idx(A::SparseMatrixCSC)
    LinearAlgebra.issymmetric(A) || throw(ArgumentError("Matrix should be symmetric"))
    n = size(A, 2)
    X = zeros(Int, nnz(A))
    rowcnt = zeros(Int, n)
    for j in 1:n
        for k in nzrange(A, j)
            r = A.rowval[k]
            rowcnt[r] += 1
            if r > j
                z = nzrange(A, r)[rowcnt[r]]
                X[k], X[z] = z, k
            end
        end
    end
    X
end

# return the index in `A.nzval` of element (`i`,`j`) if it's nonzero, 
# otherwise raise an error
#
# copied from julia/stdlib/v1.7/SparseArrays/src/sparsematrix.jl:2158
function nzindex(A::SparseMatrixCSC, i::Integer, j::Integer)
    if !(1 <= i <= size(A, 1) && 1 <= j <= size(A, 2)); throw(BoundsError()); end
    k1 = Int(A.colptr[j])
    k2 = Int(A.colptr[j+1]-1)
    (k1 > k2) && throw(ArgumentError("Matrix element ($i,$j) is zero"))  
    k = searchsortedfirst(rowvals(A), i, k1, k2, Base.Order.Forward)
    if k > k2 || A.rowval[k] != i
        throw(ArgumentError("Matrix element ($i,$j) is zero"))
    end
    return k
    # ???? Propagate errors to higher-level functions? E.g. do not throw errors
    #  here but return some Int, then from `edge_idx` print something like
    #  "Graph does not contain edge (i,j)" ?????
end

# return the indices `(i,j)` of the `k`th element in `A.nzval`
function nzindex(A::SparseMatrixCSC, k::Integer)
    j = searchsortedfirst(A.colptr, k+1) - 1
    i = rowvals(A)[k]
    i, j
end

# Void type used to fill sparse arrays w/o using memory for the nonzero values, 
#  only for their positions
struct NullNumber <: Number; end
Base.zero(::NullNumber) = false
Base.zero(::Type{NullNumber}) = false