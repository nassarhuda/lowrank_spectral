function normout(A::SparseMatrixCSC{T,Int64}) where T
  # NORMOUT Normalize the outdegrees of the matrix A.
  #
  # P = normout(A)
  #
  #   P has the same non-zero structure as A, but is normalized such that the
  #   sum of each row is 1, assuming that A has non-negative entries.
  #

  # compute the row-sums/degrees
  d = sum(A,2)
  d = squeeze(d',1)
  id = ones(size(A,1))./d
  P = spdiagm(id)*A
  return P

end
