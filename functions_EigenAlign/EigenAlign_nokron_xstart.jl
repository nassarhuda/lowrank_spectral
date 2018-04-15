# EigenAlign Algorithm
# input: two matrices A and B
# input: s1,s2,s3

#TODO: input types
#TODO: error checks

using Munkres

function EigenAlign_nokron_xstart(A,B,s1,s2,s3,iters,Xstart)
  if size(A,1)<=size(B,1)

    (ei,ej,Xmat,weight,conserved_edges) = EigenAlign_helper_nokron_xstart(A,B,s1,s2,s3,iters,Xstart)
    return (ei,ej,Xmat,weight,conserved_edges)
  else
     (ei,ej,Xmat,weight,conserved_edges) = EigenAlign_helper_nokron_xstart(B,A,s1,s2,s3,iters,Xstart)
    return (ej,ei,Xmat',weight,conserved_edges)
  end
end
function EigenAlign_helper_nokron_xstart(A,B,s1,s2,s3,iters,Xstart)
  #error checks
  gam1 = s1+s2-2s3
  gam2 = s3-s2
  gam3 = s2

  nA = size(A,1)
  nB = size(B,1)

  # AkronB = kron(A,B)
  # AkronE = kron(A,ones(nB,nB))
  # EkronB = kron(ones(nA,nA),B)
  # EkronE = kron(ones(nA,nA),ones(nB,nB))
  #
  # M = gam1*AkronB + gam2*AkronE + gam2*EkronB + gam3*EkronE

  Eb = ones(Int64,nB,nB)
  Ea = ones(Int64,nA,nA)

  # Power iteration
  # X = ones(nA*nB,1)./(nA*nB)
  # X = X./norm(X,1)

  # X = ones(nB,nA)
  # X = X./sum(X)

  X = Xstart
  X = X./sum(X)

  for i = 1:iters
    println("iteration $i started")
    # X = gam1*B*X*A' + gam2*Eb*X*A' + gam2*B*X*Ea' + gam3*Eb*X*Ea'
    X = gam1*B*X*A' + gam2*Eb*X*A' + gam2*B*X*Ea + gam3*Eb*X*Ea
    # X = M*X
    X = X./sum(X)
  end
  X
  Xmat = reshape(X,nB,nA)

  # Run Hungarian method
  # Xmat = Xmat'
  # ej = munkres(-Xmat)
  # ei = 1:length(ej)
  # ids = find(ej)
  # ej = ej[ids]
  # ei = ei[ids]

  # or bipartite matching
  Xmat = Xmat'
  ei,ej = edge_list(bipartite_matching(sparse(Xmat)))

  MATCHING = sparse(ei,ej,1,nA,nB)
  weight = X[:]'*MATCHING'[:]

  Ai = A[ei,ei]
  Bi = B[ej,ej]

  conserved_edges = nnz(Ai.*Bi)/2

  return (ei,ej,Xmat,weight,conserved_edges)

end
