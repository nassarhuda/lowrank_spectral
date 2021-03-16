# EigenAlign Algorithm
# input: two matrices A and B
# input: s1,s2,s3

# using Munkres

function EigenAlign_while_resid(A,B,s1,s2,s3,residerror)
  #error checks
  gam1 = s1+s2-2s3
  gam2 = s3-s2
  gam3 = s2

  nA = size(A,1)
  nB = size(B,1)


  Eb = ones(Int64,nB,nB)
  Ea = ones(Int64,nA,nA)


  X = ones(nB,nA)
  X = X./sum(X)

  tic()
  resid = 1
  iter = 1
  while resid >residerror
    iter+=1
    Y = gam1*B*X*A' + gam2*Eb*X*A' + gam2*B*X*Ea + gam3*Eb*X*Ea  #(Mx)
    y = Y[:]
    x = X[:]
    lam = (x'*y)./(x'*x) #(x'Mx/x'x)
    resid = norm(y - lam[1]*x)/(abs(lam[1])+1)
    X = Y./norm(Y[:],1)
  end
  toc()
  Xmat = reshape(X,nB,nA)*nB*nA;

  tic(); ej,ei = edge_list(bipartite_matching(sparse(Xmat)));t = toc()

  if resid > 1e-10
    error("residual is bigger than 1e-10")
  end

  return ei,ej,t
end
