include("readSMAT.jl")
A_dblp = readSMAT("dblp_A.smat")
A = (A_dblp'*A_dblp)
A = spones(A)
diagids = sub2ind(size(A),1:size(A,1),1:size(A,2))
A[diagids] = 0
dropzeros!(A) #avg degree 5.788682625613863

D = triu(A)
ei,ej,ev = findnz(D);


###
include("../../all_code/include_all.jl")
include(joinpath(Pkg.dir(),"handy_julia","handy_julia.jl"))
using NetworkAlign
using MatrixNetworks
using NumbersFromText
using Plots
pyplot()

Jaccard = zeros(length(ei))
for k = 1:length(ei)
  i = ei[k]
  j = ej[k]
  num = nnz(A[:,i].*A[:,j])
  den = nnz(A[:,i]+A[:,j])
  @show Jaccard[k] = num/den
end

function extract_ego(A,id)
  neighbours = A[:,id]
  loc,vals = findnz(neighbours)
  M = A[loc,loc]
  return M
end

function extract_ego_with_egonode(A,id)
  neighbours = A[:,id]
  loc,vals = findnz(neighbours)
  loc = vcat(loc,id)
  M = A[loc,loc]
  return M
end

degrees = sum(A,2)

s = sortperm(Jaccard,rev=true)
v = sort(Jaccard,rev=true)
ei_sorted = ei[s]
ej_sorted = ej[s]
scores_withego = zeros(length(ei))
scores2_withego = zeros(length(ei))
scores = zeros(length(ei))
scores2 = zeros(length(ei))
max_sizes = zeros(length(ei))
max_induced_sizes_withego = zeros(length(ei))
for i = 1:length(ei)
  @show i
  k = s[i]
  net1 = ei[k]
  net2 = ej[k]
  M = extract_ego_with_egonode(A,net1)
  N = extract_ego_with_egonode(A,net2)
  tic(); ma,mb,D = align_networks_eigenalign(M,N,8,"lowrank_svd_union",3);toc()
  Mi = M[ma,ma]
  Ni = N[mb,mb]
  c = nnz(Mi.*Ni)
  scores_withego[i] = c/max(nnz(Mi),nnz(Ni))
  scores2_withego[i] = c/max(nnz(M),nnz(N))
  max_induced_sizes_withego[i] = max(nnz(Mi),nnz(Ni))
  max_sizes[i] = max(nnz(M),nnz(N))
end

######################################
# Analysis starts here
######################################
valid_edges = zeros(Bool,length(ei))
neighbors_size = 100
for i = 1:length(ei)
  @show i
  k = s[i]
  node1 = ei[k]
  node2 = ej[k]
  if degrees[node1] >= neighbors_size && degrees[node2] >= neighbors_size
    valid_edges[i] = true
  end
end
######################################

J = Jaccard[s][valid_edges]
curscores = scores2_withego[valid_edges]
cursizes = max_induced_sizes_withego[valid_edges]

######################################

Plots.scatter(J,curscores,xlabel="Jaccard similarity score",ylabel="normalized overlap",legend=false,grid = false)
# Plots.scatter(curscores,label="normalized overlap")
# Plots.plot!(J,label="Jaccard similarity")
imgfilename = join(["dblp_jaccard_overlap",neighbors_size,"neighbors_normalized_original.pdf"])
savefig(imgfilename)

# #### Loess
# model = loess(Float64.(collect(1:length(curscores))), curscores)
# us = collect(1:0.1:length(curscores))
# vs = predict(model, us)
#
# Plots.scatter(curscores,label="normalized overlap")
# Plots.plot!(us,vs,label="Loess Fit")
# Plots.plot!(J,label="Jaccard similarity")
# imgfilename = join(["dblp_loess_fit_",neighbors_size,"_neighbors.pdf"])
# savefig(imgfilename)
