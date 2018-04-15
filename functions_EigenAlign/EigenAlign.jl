# EigenAlign Algorithm
# input: two matrices A and B
# input: s1,s2,s3

#TODO: input types
#TODO: error checks

using Munkres

function EigenAlign(A,B,s1,s2,s3,iters)
  if size(A,1)<=size(B,1)

    (ei,ej,Xmat,weight,conserved_edges) = EigenAlign_helper(A,B,s1,s2,s3,iters)
    return (ei,ej,Xmat,weight,conserved_edges)
  else
     (ei,ej,Xmat,weight,conserved_edges) = EigenAlign_helper(B,A,s1,s2,s3,iters)
    return (ej,ei,Xmat',weight,conserved_edges)
  end
end
function EigenAlign_helper(A,B,s1,s2,s3,iters)
  #error checks
  gam1 = s1+s2-2s3
  gam2 = s3-s2
  gam3 = s2

  nA = size(A,1)
  nB = size(B,1)

  AkronB = kron(A,B)
  AkronE = kron(A,ones(nB,nB))
  EkronB = kron(ones(nA,nA),B)
  EkronE = kron(ones(nA,nA),ones(nB,nB))

  M = gam1*AkronB + gam2*AkronE + gam2*EkronB + gam3*EkronE

  # Power iteration
  X = ones(nA*nB,1)./(nA*nB)
  X = X./norm(X,2)



  x = copy(X)
  for i = 1:iters
    y = M*x
    x = y./norm(y)
    lam = x'*y
    @show lam
  end
  Xmat = reshape(x,nB,nA)

  # for i = 1:iters
  #   X = M*X
  #   X = X./norm(X,2)
  # end
  # X
  # Xmat = reshape(X,nB,nA)

  # Run Hungarian method
  # Xmat = Xmat'
  # ej = munkres(-Xmat)
  # ei = 1:length(ej)
  # ids = find(ej)
  # ej = ej[ids]
  # ei = ei[ids]

  # or bipartite matching
  # try using intmatch
  Xmat = Xmat'
  ei,ej = edge_list(bipartite_matching(sparse(Xmat)))

  MATCHING = sparse(ei,ej,1,nA,nB)
  weight = X'*MATCHING'[:]

  Ai = A[ei,ei]
  Bi = B[ej,ej]

  conserved_edges = nnz(Ai.*Bi)/2

  return ei,ej,x,lam #,Xmat,weight,conserved_edges

end
