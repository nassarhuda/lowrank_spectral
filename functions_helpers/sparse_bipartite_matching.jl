function sparse_bipartite_matching(ma,mb,U,V)
  nU = size(U,1)
  nV = size(V,1)
  u = U[ma,:]
  v = V[mb,:]
  weights = vec(sum(u.*v,2))
  X = sparse(ma,mb,weights,nU,nV)
  mar,mbr = edge_list(bipartite_matching(X))
  return mar,mbr
end
