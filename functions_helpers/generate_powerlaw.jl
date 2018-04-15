# TODO error checks
# TODO types of input

function generate_powerlaw(A,n,theta)
  nA = size(A,1)
  M = spzeros(n,n)
  M[1:nA,1:nA] = A
  for i = nA+1:n
    # get degrees
    degs = sum(M[1:i-1,1:i-1],1)
    degs = degs./sum(degs)
    degs = vec(degs)
    degs = cumsum(degs)
    for j = 1:theta
      rn = rand()
      id = findin_range(rn,degs)
      M[i,id] = 1
      M[id,i] = 1
    end
  end
  return M
end

function generate_random_graph(n)
  A = sprand(n,n,0.5)
  A[1:n+1:end] = 0
  A = spones(A)
  A = convert(SparseMatrixCSC{Int64,Int64},A)
  A = max.(A,A')
  return A
end
