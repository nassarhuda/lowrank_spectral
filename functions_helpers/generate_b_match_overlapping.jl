function generate_b_match_overlapping(ma,mb,b)
  n = length(ma)
  prev = div(b,2)
  P = sprand(Int64,n,n,0.)
  for i = 1:n
    j = i-prev:i-prev+b-1
    # j = i:i+b-1
    ids = j[j.>0]
    ids = collect((ids[ids.<=n]))
    for k = ids
      P[ma[i],mb[k]] = 1
    end
  end
  return P
end

function generate_b_match_overlapping(ma,mb,b,n,nmin)
  prev = div(b,2)
  P = sprand(Int64,n,n,0.)
  for i = 1:nmin
    j = i-prev:i-prev+b-1
    # j = i:i+b-1
    ids = j[j.>0]
    ids = collect((ids[ids.<=nmin]))
    for k = ids
      P[ma[i],mb[k]] = 1
    end
  end
  return P
end
