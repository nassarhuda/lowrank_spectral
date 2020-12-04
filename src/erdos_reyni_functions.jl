using MatrixNetworks
using Random 

# this function generates a pair of erdos reyni networks to be aligned
function generate_erdos_reyni_pair(n::Int64,p::Float64,noise_model::Int64,Pe::Float64)
  #first generate A
  A = sparse(erdos_renyi_undirected(n,p))
  Q = sparse(erdos_renyi_undirected(n,Pe))
  # Q = rand(size(A))
  # Q = Q .>= Pe
  # Q = convert(SparseMatrixCSC{Int64,Int64},Q)
  if noise_model == 1
    #create B using noise model 1
    Atilde1 = A.*(1.0 .- Q) + (1.0 .- A).*Q # noise model 1
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
    Atilde2 = A.*(1.0 .-Q) + (1.0 .-A).*Qtil # noise model 2
    P = sparse(collect(1:n),collect(n:-1:1),1)
    B = P*Atilde2*P'
  elseif noise_model == 3
    #create B using noise model 1
    Atilde1 = A.*(1.0 .- Q) + (1.0 .- A).*Q # noise model 1
    # P = sparse(randperm(nA),randperm(nA),1)
    perm = shuffle(1:n)
    P = sparse(collect(1:n),perm,1)
    B = P*Atilde1*P'
    return A,B,perm
  else
    error("noise_model should either be 1 or 2")
  end
  return A,B
end

# methods:
# netalignmr
include("/Users/ccolley/Code/NetworkAlignment.jl/src/NetworkAlignment.jl") #local NetworkAlignment code
function align_erods_reyni(A,B,method::String)
  L = sparse(ones(Int64,size(A,1),size(B,1)))
  S,w,li,lj = netalign_setup(A,B,L)

  # align networks
  a = 1;
  b = 1;
  stepm = 25;
  rtype = 1;
  maxiter = 10;
  verbose = true;
  gamma = 0.4;

  if method == "netalignmr"
    xbest,st,status,hist = netalignmr(S,w,a,b,li,lj,gamma,stepm,rtype,maxiter,verbose)
  elseif method == "isorank"
    xbest,flag,reshist = isorank(S,w,a,b,li,lj,b/(a+b),2,1e-12,maxiter)
  elseif method == "netalignbp"
    xbest,flag,reshist = netalignbp(S,w,a,b,li,lj,gamma,1,maxiter,verbose)
    # x,flag,reshist = netalignbp(S,w,a,b,li,lj,0.99,2,100,true)
  else
    error("Alignment method should be either (1) \"netalignmr\" or (2) \"netalignbp\" or (3) \"isorank\"")
  end

  BM = bipartite_matching(xbest,li,lj)
  ma,mb = edge_list(BM)
  return ma,mb
end
# evaluation:
function evaluate_erdosreyni_experiment(A,B,ma,mb,perm=nothing)
  nA = size(A,1)

  if perm === nothing
    P = sparse(collect(1:nA),collect(nA:-1:1),1)
  else
    P = sparse(collect(1:nA),perm,1)
  end
  Ptil = sparse(ma,mb,1,nA,nA)
  recov = 1-(norm(P-Ptil,1)/(2nA))
  return recov
end
