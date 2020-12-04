function generate_b_match_overlapping(ma,mb,b,nU,nV)
  n = length(mb)
  ln = length(ma)
  prev = div(b,2)
  P = spzeros(Int64,nU,nV)
  for i = 1:ln
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
