function newbound_rounding(U,V,bmatchval)
  U_sortperm = sortcolsperm(U,true)
  V_sortperm = sortcolsperm(V,true)
  nU = size(U,1)
  nV = size(V,1)
  r = size(U,2)
  @assert r==size(V,2)

  d = min(size(U_sortperm,1),size(V_sortperm,1))
  U_weights = sort(U,1,rev=true)
  V_weights = sort(V,1,rev=true)

#   U_weights = U_weights[1:d,:]
#   V_weights = V_weights[1:d,:]
#
#   U_sortperm = U_sortperm[1:d,:]
#   V_sortperm = V_sortperm[1:d,:]


  P = spzeros(nU,nV)
  allrecoveries = zeros(size(U,2))
  for i = 1:size(U,2)
    ui = U_weights[:,i]
    vi = V_weights[:,i]
    lastid_ui = findfirst(ui.<0)
    lastid_vi = findfirst(vi.<0)

    if lastid_ui == 0 && lastid_vi == 0
      lastidpos = d
      lneg = -1
    elseif lastid_vi == 0
      lastidpos = min(d,lastid_ui-1)
      lneg = -1
    elseif lastid_ui == 0
      lastidpos = min(d,lastid_vi-1)
      lneg = -1
    else
      lastidpos = min(lastid_ui,lastid_vi)-1
      lneg = min(nU-lastid_ui,nV-lastid_vi)
    end

    ei1 = U_sortperm[1:lastidpos,i]
    ej1 = V_sortperm[1:lastidpos,i]
    ei2 = U_sortperm[nU-lneg:nU,i]
    ej2 = V_sortperm[nV-lneg:nV,i]
    ei = vcat(ei1,ei2)
    ej = vcat(ej1,ej2)
    # allrecoveries[i] = evaluate_erdosreyni_experiment(A,B,ei,ej)
    P = P + sparse(ei,ej,1,nU,nV) #+ generate_b_match_overlapping(ei,ej,2,nU,nV)
    # with bmatching:
    # bmatchval = 20
    for bm = 1:bmatchval
    if !isempty(ej)
      popfirst!(ej)
      P = P + sparse(ei[1:end-bm],ej,1,nU,nV)
    end
    end


  end
  @show nnz(P)/prod(size(P))
  return P
end

function newbound_rounding_lowrank_evaluation_relaxed(U,V,bmatchval)
  # U = Float32.(U)
  # V = Float32.(V)
  U_sortperm = sortcolsperm(U,true)
  V_sortperm = sortcolsperm(V,true)
  nU = size(U,1)
  nV = size(V,1)
  r = size(U,2)
  @assert r == size(V,2)

  d = min(size(U_sortperm,1),size(V_sortperm,1))
  U_weights = sort(U,dims=1,rev=true)
  V_weights = sort(V,dims=1,rev=true)

#   U_weights = U_weights[1:d,:]
#   V_weights = V_weights[1:d,:]
#
#   U_sortperm = U_sortperm[1:d,:]
#   V_sortperm = V_sortperm[1:d,:]


#   P = spzeros(nU,nV)
  U1 = []
  V1 = []
  allrecoveries = zeros(size(U,2))
  for i = 1:size(U,2)
    ui = U_weights[:,i]
    vi = V_weights[:,i]
    lastid_ui = findfirst(ui.<0)
    lastid_vi = findfirst(vi.<0)

    lastid_ui = (lastid_ui === nothing) ? 0 : lastid_ui
    lastid_vi = (lastid_vi === nothing) ? 0 : lastid_vi

    if lastid_ui == 0 && lastid_vi == 0
      lastidpos = d
      lneg = -1
    elseif lastid_vi == 0
      lastidpos = min(d,lastid_ui-1)
      lneg = -1
    elseif lastid_ui == 0
      lastidpos = min(d,lastid_vi-1)
      lneg = -1
    else
      lastidpos = min(lastid_ui,lastid_vi)-1
      lneg = min(nU-lastid_ui,nV-lastid_vi)
    end

    ei1 = U_sortperm[1:lastidpos,i]
    ej1 = V_sortperm[1:lastidpos,i]
    ei2 = U_sortperm[nU-lneg:nU,i]
    ej2 = V_sortperm[nV-lneg:nV,i]
    ei = vcat(ei1,ei2)
    ej = vcat(ej1,ej2)
    # allrecoveries[i] = evaluate_erdosreyni_experiment(A,B,ei,ej)
#     P = P + sparse(ei,ej,1,nU,nV)

    # P = P + sparse(ei,ej,1,nU,nV) #+ generate_b_match_overlapping(ei,ej,2,nU,nV)
    # with bmatching:
    append!(U1,ei)
    append!(V1,ej)
    # bmatchval = 6
    #println("bmatchval is $bmatchval")
    for bm = 1:bmatchval
      if !isempty(ej)
        popfirst!(ej)
        # P = P + sparse(ei[1:end-bm],ej,1,nU,nV)
        append!(U1,ei[1:end-bm])
        append!(V1,ej)
      end
    end
  end

  all_matches = [U1 V1]
  unique_matches = unique(all_matches,dims=1)
  U1unique = unique_matches[:,1]
  V1unique = unique_matches[:,2]
  uo = U[U1unique,:]
  vo = V[V1unique,:]
  weights = vec(sum(uo.*vo,dims=2))
  X = sparse(U1unique,V1unique,weights,nU,nV)

  return X
end
