function read_net_file(filename::String)
  M = readdlm(filename)
  ei = M[:,1]
  ej = M[:,2]

  ei = map(x->x[2:end],ei)
  ej = map(x->x[2:end],ej)

  ei = map(x->parse(Int64,x),ei)
  ej = map(x->parse(Int64,x),ej)

  m = maximum(ei)
  n = maximum(ej)
  d = max(m,n)

  return sparse(ei,ej,1,d,d)
end

function read_sim_file(filename::String,m::Int64,n::Int64)
  M = readdlm(filename)
  ei = M[:,1]
  ej = M[:,2]

  ei = map(x->x[2:end],ei)
  ej = map(x->x[2:end],ej)

  ei = map(x->parse(Int64,x),ei)
  ej = map(x->parse(Int64,x),ej)

  w = M[:,3]
  w = convert(Array{Float64},w)

  return sparse(ei,ej,w,m,n)
end

function setup_yeast_networks(file1,file2)
    lb1,A1 = convert_LEDA_to_adjacency(file1)
    lb2,A2 = convert_LEDA_to_adjacency(file2)
    n = size(A2,1)
    ei,ej,ev = findnz(A2)
    eii = n+1-ei # flip to unbias
    ejj = n+1-ej # flip to unbias
    B = sparse(eii,ejj,ev)
    Pref = sparse(1:n,collect(n:-1:1),1,size(A1,1),size(B,1))
    return A1,B,Pref
end

function gs3_evaluate(A,B,mi,mj)
  mi = collect(mi)
  mj = collect(mj)
  num_matches = length(mi)
  NCV = (length(unique(mi)) + length(unique(mj))) / (size(A, 1) + size(B, 1));

  #ii,jj,vv = findnz(A)
  #ei = mi[ii]
  #ej = mi[jj]
  #A_induced = sparse(ei,ej,vv)

  #ii,jj,vv = findnz(B)
  #ei = mj[ii]
  #ej = mj[jj]
  #B_induced = sparse(ei,ej,vv)

  A_induced = A[mi, mi];
  B_induced = B[mj, mj];
  D = A_induced + B_induced;
  D = D - spdiagm(diag(D))
  ii,jj = ind2sub(D,find(D.==2))
  align_graph = sparse(ii, jj, 1, num_matches, num_matches)

  M_align = nnz(align_graph) / 2
  # T_align = full(trace(align_graph^3)/6);

  GS3 = M_align / (nnz(A_induced)/2 + nnz(B_induced)/2 - M_align)
  NCV_GS3 = sqrt(GS3*NCV)
  return NCV_GS3
end

function gs3_evaluate2(A,B,mi,mj)
  mi = collect(mi)
  mj = collect(mj)
  num_matches = length(mi)
  NCV = (length(unique(mi)) + length(unique(mj))) / (size(A, 1) + size(B, 1));

  #ii,jj,vv = findnz(A)
  #ei = mi[ii]
  #ej = mi[jj]
  #A_induced = sparse(ei,ej,vv)

  #ii,jj,vv = findnz(B)
  #ei = mj[ii]
  #ej = mj[jj]
  #B_induced = sparse(ei,ej,vv)

  A_induced = A[mi, mi];
  B_induced = B[mj, mj];
  D = A_induced + B_induced;
  D = D - spdiagm(diag(D))
  ii,jj = ind2sub(D,find(D.==2))
  align_graph = sparse(ii, jj, 1, num_matches, num_matches)

  M_align = nnz(align_graph) / 2
  # T_align = full(trace(align_graph^3)/6);

  GS3 = M_align / (nnz(A_induced)/2 + nnz(B_induced)/2 - M_align)
  NCV_GS3 = 1/(((1/GS3)+(1/NCV))/2)
  return NCV_GS3
end

"""
sortcolsperm(A,true) returns the indices of sorted columns in A in descending order
sortcolsperm(A,false) returns the indices of sorted columns in A in ascending order
```
example:
  julia> W = rand(3,3)
    3×3 Array{Float64,2}:
    0.661943  0.00517749  0.332394
    0.716344  0.61179     0.544258
    0.372336  0.994069    0.297704
  julia> sortcolsperm(W,true)
    3×3 Array{Int64,2}:
    2  3  2
    1  2  1
    3  1  3
  julia> sortcolsperm(W,false)
    3×3 Array{Int64,2}:
    3  1  3
    1  2  1
    2  3  2
```
"""
function sortcolsperm(X::Matrix{T},REV::Bool) where T
    P = Matrix{Int}(undef,size(X,1),size(X,2))

    Threads.@threads for i=1:size(X,2)
        P[:,i] = sortperm(X[:,i]; rev=REV)
    end
    return P
end
