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
  lastid = n;

  for i = 1:n
    curmax = 0
    curma = 0
    curmb = 0
    for mi = ma_ids
      u = U1t[:,mi]
      for mj = mb_ids

        v = V1t[:,mj]
        m = dot(u,v)
        if m > curmax
          curmax = m
          curma = mi
          curmb = mj
        end
      end
    end

    if curma == 0
      lastid = i-1
      break;
    end
    @show curma
    @show curmb

    ma[i] = curma
    mb[i] = curmb
    ma_ids = setdiff(ma_ids,curma)
    mb_ids = setdiff(mb_ids,curmb)
  end
  ma = ma[1:lastid]
  mb = mb[1:lastid]

  return ma,mb
end

function greedy_lowrank(U1,V1)

  k = size(U1,2)
  n = size(U1,1)

  ma = zeros(Int64,n)
  mb = zeros(Int64,n)

  ma_ids = collect(1:n)
  mb_ids = collect(1:n)

  U1t = U1';
  V1t = V1';
  lastid = n;

  for i = 1:n
    curmax = 0
    curma = 0
    curmb = 0
    for mi = ma_ids
      u = U1t[:,mi]
      for mj = mb_ids

        v = V1t[:,mj]
        m = dot(u,v)
        if m > curmax
          curmax = m
          curma = mi
          curmb = mj
        end
      end
    end

    if curma == 0
      lastid = i-1
      break;
    end
    @show curma
    @show curmb
    ma[i] = curma
    mb[i] = curmb
    ma_ids = setdiff(ma_ids,curma)
    mb_ids = setdiff(mb_ids,curmb)
  end

  ma = ma[1:lastid]
  mb = mb[1:lastid]

  return ma,mb
end
