# sparsified greedy matching
include("decomposeX.jl")
using MatrixNetworks

function keep_top_k(v,topk)
  vsortedperm = sortperm(v,rev=true)
  tokeep_ids = vsortedperm[1:topk]
  tokeep_vals = v[tokeep_ids]
  n = length(v)
  w = sprand(n,1,0.0)
  w[tokeep_ids] = tokeep_vals
  return w
end

function sparsified_EigenAlign(A,B,c1,c2,c3,iters,topk)
  Uk,Vk,Wk,W1,W2 = decomposeX(A,B,iters,c1,c2,c3)
  U1 = Uk
  V1 = Vk*Wk'

  k = size(U1,2)
  n = size(U1,1)

  X = spzeros(n,n)
  for rowid = 1:n
    u = U1[rowid,:]
    tempvec = zeros(n)
    for i = 1:length(u)
      tempvec += u[i] * V1[:,i]
    end
    w = keep_top_k(tempvec,topk)
    X[rowid,find(w)] = w[find(w)]
  end
  ma,mb = edge_list(bipartite_matching(X))

  return ma,mb,X
end
