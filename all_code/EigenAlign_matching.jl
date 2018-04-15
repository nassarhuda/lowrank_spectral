# EigenAlign Algorithm
# input: two matrices A and B
# input: s1,s2,s3

#TODO: input types
#TODO: error checks

using Munkres

function EigenAlign_matching(A,B,s1,s2,s3,iters)
  if size(A,1)<=size(B,1)

    ei,ej = EigenAlign_helper_matching(A,B,s1,s2,s3,iters)
    return ei,ej
  else
     ei,ej = EigenAlign_helper_matching(B,A,s1,s2,s3,iters)
    return ej,ei
  end
end
function EigenAlign_helper_matching(A,B,s1,s2,s3,iters)

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

  X = ones(nB,nA)
  X = X./sum(X)

  for i = 1:iters
    X = gam1*B*X*A' + gam2*Eb*X*A' + gam2*B*X*Ea + gam3*Eb*X*Ea
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
  ej,ei = edge_list(bipartite_matching(sparse(Xmat)))
  ## bmatching start
  P = sparse(ei,ej,1,nA,nB)
  shift!(ej)
  P = P + sparse(ei[1:end-1],ej,1,nA,nB)
  Xsample = P'.*Xmat
  ej,ei = edge_list(bipartite_matching(sparse(Xsample)))
  ## bmatching end

  return ei,ej

end
