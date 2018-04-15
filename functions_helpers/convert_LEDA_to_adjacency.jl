function convert_LEDA_to_adjacency(filename::String)
  M = readdlm(filename)
  nnodes = M[4,1]
  nodelabels = M[4+(1:nnodes),1]

  nedges = M[4+nnodes+1,1]
  edges = M[4+nnodes+1+(1:nedges),1:2]

  A = sparse(edges[:,1],edges[:,2],1,nnodes,nnodes)
  A = max.(A,A')

  return nodelabels,A
end
