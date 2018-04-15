function generate_b_match(ma,mb,b)
  n = max(maximum(ma),maximum(mb))
  P = zeros(Int64,n,n)
  ln = length(ma)
  remainder = mod(ln,b)
  ls = ln-remainder
  for i = 1:b:ls

    ids = i:i+b-1
    P[ma[ids],mb[ids]] = 1
  end
  ids = ls+1:n
  P[ids,ids] = 1
  return P
end
