using MatrixNetworks

function generate_powerlaw(n,theta)
  assert(n>=5)
  M = spzeros(Int64,n,n)
  if theta == 4
    netden = 0.08
  elseif theta == 6
    netden = 0.11
  end
  A = sparse(erdos_renyi_undirected(6,netden))
  while nnz(A) == 0
    A = sparse(erdos_renyi_undirected(6,netden))
  end
  # A = rand(6,6)
  # A = A .>= 0.5
  # A = A.*A'
  A = convert(SparseMatrixCSC{Int64,Int64},A)
  M[1:6,1:6] = A
  for i = 7:n
    A = M[1:i-1,1:i-1]
    deg = sum(A,2)
    degp = deg./sum(deg)
    c = cumsum(degp)
    r = rand(theta)
    nodes = map(x->findfirst(c.>x),r)
    M[nodes,i] = 1
    M[i,nodes] = 1
  end
  return M
end

# this function generates a pair of erdos reyni networks to be aligned
function generate_powerlaw_pair(n::Int64,theta::Int64,noise_model::Int64,Pe::Float64)
  #first generate A
  A = generate_powerlaw(n,theta)
  p = nnz(A)/prod(size(A))
  Q = sparse(erdos_renyi_undirected(n,Pe))
  # Q = rand(size(A))
  # Q = Q .>= Pe
  # Q = convert(SparseMatrixCSC{Int64,Int64},Q)
  if noise_model == 1
    #create B using noise model 1
    Atilde1 = A.*(1-Q) + (1-A).*Q # noise model 1
    # P = sparse(randperm(nA),randperm(nA),1)
    P = sparse(collect(1:n),collect(n:-1:1),1)
    B = P*Atilde1*P'
  elseif noise_model == 2
    #create B using noise model 2
    Pe2 = (p*Pe)/(1-p)
    Qtil = sparse(erdos_renyi_undirected(n,Pe2))
    # Qtil = rand(size(A))
    # Qtil[Qtil.<Pe2] = 1
    # Qtil[Qtil.>=Pe2] = 0
    # Qtil = convert(SparseMatrixCSC{Int64,Int64},Qtil)
    Atilde2 = A.*(1-Q) + (1-A).*Qtil # noise model 2
    P = sparse(collect(1:n),collect(n:-1:1),1)
    B = P*Atilde2*P'
  else
    error("noise_model should either be 1 or 2")
  end
  return A,B
end
