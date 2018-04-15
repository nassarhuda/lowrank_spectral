using MatrixNetworks
include("../functions_decompositions/decomposeX.jl")
include("newbound_methods.jl")
include("greedy_lowrank_EigenAlign.jl")

function align_networks_eigenalign(A,B,iters,method,bmatch; default_params = true)
  D = 0
  s1,s2,s3 = find_parameters(A,B)
  if default_params == false
    s1 += 100
    s2 += 10
    s3 += 5
    @show s1
    @show s2
    @show s3
  end

  c1 = s1+s2-2s3; c2 = s3-s2; c3 = s2
  Uk,Vk,Wk,W1,W2 = decomposeX_balance_allfactors(A,B,iters+1,c1,c2,c3)
  Un,Vn = split_balanced_decomposition(Uk,Wk,Vk)
  timematching = 0

  nA = size(A,1)
  nB = size(B,1)

  if method == "full_multiplication"
    Xsample = 1e6*Uk*Wk*Vk'
    ma,mb = edge_list(bipartite_matching(sparse(Xsample)))
  elseif method == "lowrank_unbalanced_union"
    U = copy(Uk)
    V = Vk*Wk'
    P = newbound_rounding(U,V)
    Xsample = spones(P).*Xn
    ma,mb = edge_list(bipartite_matching(Xsample))
  elseif method == "lowrank_unbalanced_best"
    U = copy(Uk)
    V = Vk*Wk'
    ma,mb,D = newbound_rounding_max_approx(U,V)
  elseif method == "lowrank_balanced_union"
    U = Un
    V = Vn
    Xn = Un*Vn'
    P = newbound_rounding(U,V,bmatch)
    Xsample = spones(P).*Xn
    ma,mb = edge_list(bipartite_matching(Xsample))
    D = nnz(P)/prod(size(P))
  elseif method == "lowrank_balanced_union_lowrank_evaluation"
    U = Un
    V = Vn
    X = newbound_rounding_lowrank_evaluation_relaxed(U,V,bmatch).*(10^8) #bmatch
    D = nnz(X)/prod(size(X))
    ma,mb = edge_list(bipartite_matching(X))
  elseif method == "lowrank_balanced_best"
    U = Un
    V = Vn
    ma,mb,D = newbound_rounding_max_approx(U,V)
  elseif method == "lowrank_balanced_best_eff"
    U = Un
    V = Vn
    ma,mb,dval = newbound_rounding_max_approx_eff(U,V)
    X = U*V'
    X = X*100
    ei,ej = edge_list(bipartite_matching(sparse(X)))
    v0 = sum(sparse(ei,ej,1,size(X,1),size(X,2)).*X)
    M = sparse(ma,mb,1,size(X,1),size(X,2))
    v1 = sum(M.*X)
    row,col,M = greedy_match(X)
    v2 = sum(M.*X)
    D = [dval,v0/v1,v0/v2]
  elseif method == "lowrank_Wkdecomposed_union"
    U = Uk*W1
    V = Vk*W2'
    P = newbound_rounding(U,V)
    Xn = Un*Vn'
    Xsample = spones(P).*Xn
    ma,mb = edge_list(bipartite_matching(Xsample))
  elseif method == "lowrank_Wkdecomposed_best"
    U = Uk*W1
    V = Vk*W2'
    ma,mb,D = newbound_rounding_max_approx(U,V)
  elseif method == "lowrank_svd_union"

    U,S,V = svd(Wk)
    U1 = Uk*U*diagm(sqrt(S)); V1 = Vk*V*diagm(sqrt(S));
    X = newbound_rounding_lowrank_evaluation_relaxed(U1,V1,bmatch).*(10^8); #bmatch
    avgdeg = map(x->sum(X[x,:].!=0),1:size(X,1))
    avgdeg = mean(avgdeg)
    tic(); ma,mb = edge_list(bipartite_matching(X));timematching = toq()
    D = avgdeg;#nnz(X)/prod(size(X))
  elseif method == "lowrank_lu_union"
    L,U = lu(Wk,Val{false})
    U1 = Uk*L;
    V1 = Vk*U';
    X = newbound_rounding_lowrank_evaluation_relaxed(U1,V1,bmatch).*(10^8) #bmatch
    ma,mb = edge_list(bipartite_matching(X))

  elseif method == "greedy"
    U,S,V = svd(Wk)
    U1 = Uk*U*diagm(sqrt(S)); V1 = Vk*V*diagm(sqrt(S));
    ma,mb = greedy_lowrank(U1,V1)
  else
    error("method should be one of the following: (1)eigenalign,
          (2)lowrank_unbalanced_best, (3)lowrank_unbalanced_union,
          (4)lowrank_balanced_best, (5)lowrank_balanced_union,
          (6)lowrank_Wkdecomposed_best, (7)lowrank_Wkdecomposed_union")
  end

  return ma,mb,D,timematching
end

function find_parameters(A,B)
  nB = size(B,1)
  nA = size(A,1)
  nmatches = sum(A)*sum(B)
  nmismatches = sum(A)*(nB^2 - sum(B)) + sum(B)*(nA^2 - sum(A))
  mygamma = nmatches/nmismatches
  myalpha = (1/mygamma) + 1
  myeps = 0.001
  s1 = myalpha + myeps
  s2 = 1+myeps
  s3 = myeps
  return s1,s2,s3
end

function balance_decomposition(Uk,Wk,Vk)
  Du = diagm(vec(maximum(Uk,1)))
  Dv = diagm(vec(maximum(Vk,1)))
  Ukt = Uk*diagm(vec(1./maximum(Uk,1)))
  Vkt = Vk*diagm(vec(1./maximum(Vk,1)))
  rho = (maximum(Du)*maximum(Dv))
  Wkt2 = Du*Wk*Dv./rho
  rho *= maximum(Wkt2)
  Wkt2 = Wkt2/maximum(Wkt2)
  L,U = lu(Wkt2,Val{false})
  Ud = diagm(sqrt.(abs.(diag(U)))) ##here
  L2 = L*Ud
  U2 = diagm(1./sqrt.(diag(U)))*U
  Un = Ukt*L2
  Vn = Vkt*U2'
  Xn = Ukt*L2*U2*Vkt'
  return Un,Vn,Xn
end

# can use sqrt here (like above)
function split_balanced_decomposition(Uk,Wk,Vk)
  L,U = lu(Wk,Val{false})
  Ud = Diagonal(sqrt.(abs.(diag(U))))
  L2 = L*Ud
  U2 = Diagonal(1./sqrt.(diag(U)))*U
  Un = Uk*L2
  Vn = Vk*U2'
  return Un,Vn
end

function split_svd(Uk,Wk,Vk)
  U,S,V = svd(Wk)
  D = Diagonal(sqrt(S))
  Unew = Uk*U*D
  Vnew = Vk*V*D
  return Unew,Vnew
end
