# EigenAlign Algorithm
# input: two matrices A and B
# input: s1,s2,s3

#TODO: input types
#TODO: error checks

using Munkres

function EigenAlign_nokron(A,B,s1,s2,s3,iters)
  if size(A,1)<=size(B,1)

    (ei,ej,Xmat,weight,conserved_edges) = EigenAlign_helper_nokron(A,B,s1,s2,s3,iters)
    return (ei,ej,Xmat,weight,conserved_edges)
  else
     (ei,ej,Xmat,weight,conserved_edges) = EigenAlign_helper_nokron(B,A,s1,s2,s3,iters)
    return (ej,ei,Xmat',weight,conserved_edges)
  end
end
function EigenAlign_helper_nokron(A,B,s1,s2,s3,iters)
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

  X = ones(nB,nA)
  X = X./sum(X)

  # x = copy(X)
  # for i = 1:iters
  #   y = M*x
  #   x = y./norm(y)
  #   lam = x'*y
  #   @show lam
  # end
  # Xmat = reshape(x,nB,nA)
  resid = 1
  for i = 1:iters
  i = 1
  # while resid >1e-14
    @show i
    i += 1
    # println("iteration $i started")
    # X = gam1*B*X*A' + gam2*Eb*X*A' + gam2*B*X*Ea' + gam3*Eb*X*Ea'
    Y = gam1*B*X*A' + gam2*Eb*X*A' + gam2*B*X*Ea + gam3*Eb*X*Ea
    X = Y./norm(Y[:],1)

    # Y = gam1*B*X*A' + gam2*Eb*X*A' + gam2*B*X*Ea + gam3*Eb*X*Ea #(Mx)
    y = Y[:]
    x = X[:]
    lam = (x'*y)./(x'*x)
    resid = norm(y - lam[1]*x)
    # if resid < 1e-16
    #   @show resid
    # end

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
  # Xmat = Xmat'
  # Xmat = Xmat'
  # ej,ei = edge_list(bipartite_matching(sparse(Xmat*100)))

  ej,ei,M = greedy_match(Xmat)
  ejbp,eibp = edge_list(bipartite_matching(sparse(Xmat*10^(abs(log10(maximum(Xmat)))))))

  MATCHING = sparse(ei,ej,1,nA,nB)
  weight = X[:]'*MATCHING'[:]

  Ai = A[ei,ei]
  Bi = B[ej,ej]

  conserved_edges = nnz(Ai.*Bi)/2

  if resid > 1e-10
    error("residual is bigger than 1e-14")
  end

  return (ei,ej,Xmat,eibp,ejbp)

end
