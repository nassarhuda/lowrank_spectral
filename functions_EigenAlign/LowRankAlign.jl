function LowRankAlign(A,B,k)

  Aeigs = eigs(A;nev = k)
  Beigs = eigs(B;nev = k)

  V = Aeigs[2]
  U = Beigs[2]

  # can do a better job with using the eigen values instead?
  M = A*V*U'*B
  ei,ej = edge_list(bipartite_matching(sparse(M)))
  return ei,ej

end
