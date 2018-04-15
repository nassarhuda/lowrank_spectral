include("../functions_decompositions/decomposeX.jl")

function greedy_lowrank_EigenAlign(A,B,c1,c2,c3,iters)
  Uk,Vk,Wk,W1,W2 = decomposeX_balance_allfactors(A,B,iters,c1,c2,c3)
  U1,V1 = split_balanced_decomposition(Uk,Wk,Vk);

  k = size(U1,2)
  n = size(U1,1)

  ma = zeros(Int64,n)
  mb = zeros(Int64,n)

  ma_ids = collect(1:n)
  mb_ids = collect(1:n)

  U1t = U1';
  V1t = V1';

  for i = 1:n
    curmax = 0
    curma = 0
    curmb = 0
    for mi = 1:length(ma_ids)
      u = U1t[:,ma_ids[mi]]
      for mj = 1:length(mb_ids)
        v = V1t[:,mb_ids[mj]]
        m = dot(u,v)
        if m > curmax
          curmax = m
          curma = mi
          curmb = mj
        end
      end
    end
    ma[i] = ma_ids[curma]
    mb[i] = mb_ids[curmb]
    deleteat!(ma_ids, curma)
    deleteat!(mb_ids, curmb)
  end

  return ma,mb
end
