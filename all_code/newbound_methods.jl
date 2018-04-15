function newbound_rounding_details(U,V,Xn)
  U_sortperm = sortcolsperm(U,true)
  V_sortperm = sortcolsperm(V,true)
  nU = size(U,1)
  nV = size(V,1)
  r = size(U,2)
  assert(r==size(V,2))

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
  allrecoveries_accum = zeros(size(U,2))
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
    allrecoveries[i] = evaluate_erdosreyni_experiment(A,B,ei,ej)
    P = P + sparse(ei,ej,1,n,n)

    Xsample = spones(P).*Xn
    maunion,mbunion = edge_list(bipartite_matching(Xsample))
    allrecoveries_accum[i] = evaluate_erdosreyni_experiment(A,B,maunion,mbunion)

  end
  return P,allrecoveries_accum,allrecoveries
end

function newbound_rounding_max_approx_eff(U,V)
  nU = size(U,1)
  nV = size(V,1)
  nU = size(U,1)
  nV = size(V,1)
  U_sortperm = sortcolsperm(U,true)
  V_sortperm = sortcolsperm(V,true)

  r = size(U,2)
  assert(r==size(V,2))

  d = min(size(U_sortperm,1),size(V_sortperm,1))
  U_weights = sort(U,1,rev=true)
  V_weights = sort(V,1,rev=true)

#   U_weights = U_weights[1:d,:]
#   V_weights = V_weights[1:d,:]
#
#   U_sortperm = U_sortperm[1:d,:]
#   V_sortperm = V_sortperm[1:d,:]

  weights = zeros(size(U,2))
  di = zeros(Float64,size(U,2))
  Umatching = Array{Int64,1}[]
  Vmatching = Array{Int64,1}[]

  for i = 1:size(U,2)
    U1 = Int64[]
    V1 = Int64[]
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

    wi = U_weights[1:lastidpos,i]
    wj = V_weights[1:lastidpos,i]
    weights[i] += sum(wi.*wj)
    wi = U_weights[nU-lneg:nU,i]
    wj = V_weights[nV-lneg:nV,i]
    weights[i] += sum(wi.*wj)

    ei = U_sortperm[1:lastidpos,i]
    ej = V_sortperm[1:lastidpos,i]
    append!(U1,ei)
    append!(V1,ej)
    ei = U_sortperm[nU-lneg:nU,i]
    ej = V_sortperm[nV-lneg:nV,i]
    append!(U1,ei)
    append!(V1,ej)

    push!(Umatching,U1)
    push!(Vmatching,V1)
  end


  for i = 1:size(U,2)
    ei = Umatching[i]
    ej = Vmatching[i]

    for j = 1:size(U,2)
      num = weights[j]
      ui = U[:,j]
      vi = V[:,j]
      ui = ui[ei]
      vi = vi[ej]
      den = sum(ui.*vi)
      f = num/den
      di[i] = max(di[i],f)
    end
  end


  D = minimum(di)
  Dind = indmin(di)

  utouse = Umatching[Dind]
  vtouse = Vmatching[Dind]

  return utouse,vtouse,D
end

function newbound_rounding_max_approx(U,V)
  nU = size(U,1)
  nV = size(V,1)
  n = size(U,1)
  U_sortperm = sortcolsperm(U,true)
  V_sortperm = sortcolsperm(V,true)

  r = size(U,2)
  assert(r==size(V,2))

  d = min(size(U_sortperm,1),size(V_sortperm,1))
  U_weights = sort(U,1,rev=true)
  V_weights = sort(V,1,rev=true)

#   U_weights = U_weights[1:d,:]
#   V_weights = V_weights[1:d,:]
#
#   U_sortperm = U_sortperm[1:d,:]
#   V_sortperm = V_sortperm[1:d,:]

  weights = zeros(size(U,2))
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

    wi = U_weights[1:lastidpos,i]
    wj = V_weights[1:lastidpos,i]
    weights[i] += sum(wi.*wj)
    wi = U_weights[nU-lneg:nU,i]
    wj = V_weights[nV-lneg:nV,i]
    weights[i] += sum(wi.*wj)
  end

  di = zeros(Float64,size(U,2))
  for i = 1:size(U,2)

    P = spzeros(nU,nV)
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

    ei = U_sortperm[1:lastidpos,i]
    ej = V_sortperm[1:lastidpos,i]
    P = P + sparse(ei,ej,1,nU,nV)
    ei = U_sortperm[nU-lneg:nU,i]
    ej = V_sortperm[nV-lneg:nV,i]
    P = P + sparse(ei,ej,1,nU,nV)

    for j = 1:size(U,2)
      num = weights[j]
      ui = U[:,j]
      vi = V[:,j]
      M = ui*vi'
      den = sum(spones(P).*M)
      f = num/den
      di[i] = max(di[i],f)
    end
  end

  D = minimum(di)
  Dind = indmin(di)

  ui = U_weights[:,Dind]
  vi = V_weights[:,Dind]

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

  ei = U_sortperm[1:lastidpos,Dind]
  ej = V_sortperm[1:lastidpos,Dind]
  ei2 = U_sortperm[nU-lneg:nU,Dind]
  ej2 = V_sortperm[nV-lneg:nV,Dind]

  utouse = vcat(ei,ei2)
  vtouse = vcat(ej,ej2)

  return utouse,vtouse,D
end


function newbound_rounding(U,V)
  U_sortperm = sortcolsperm(U,true)
  V_sortperm = sortcolsperm(V,true)
  nU = size(U,1)
  nV = size(V,1)
  r = size(U,2)
  assert(r==size(V,2))

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

    # with bmatching: moved to newbound_rounding(U,V,bmatchval)
    # bmatchval = 20
    # for bm = 1:bmatchval
    # if !isempty(ej)
    #   shift!(ej)
    #   P = P + sparse(ei[1:end-bm],ej,1,nU,nV)
    # end
    # end


  end
  @show nnz(P)/prod(size(P))
  return P
end

function newbound_rounding_lowrank_evaluation(U,V)
  U_sortperm = sortcolsperm(U,true)
  V_sortperm = sortcolsperm(V,true)
  nU = size(U,1)
  nV = size(V,1)
  r = size(U,2)
  assert(r==size(V,2))

  d = min(size(U_sortperm,1),size(V_sortperm,1))
  U_weights = sort(U,1,rev=true)
  V_weights = sort(V,1,rev=true)

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
    append!(U1,ei)
    append!(V1,ej)
  end

  all_matches = [U1 V1]
  unique_matches = unique(all_matches,1)
  U1unique = unique_matches[:,1]
  V1unique = unique_matches[:,2]
  uo = U[U1unique,:]
  vo = V[V1unique,:]
  weights = vec(sum(uo.*vo,2))
  X = sparse(U1unique,V1unique,weights,nU,nV)

  return X
end

function newbound_rounding_lowrank_evaluation_relaxed(U,V)
  U_sortperm = sortcolsperm(U,true)
  V_sortperm = sortcolsperm(V,true)
  nU = size(U,1)
  nV = size(V,1)
  r = size(U,2)
  assert(r==size(V,2))

  d = min(size(U_sortperm,1),size(V_sortperm,1))
  U_weights = sort(U,1,rev=true)
  V_weights = sort(V,1,rev=true)

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
    bmatchval = 6
    println("bmatchval is $bmatchval")
    for bm = 1:bmatchval
      if !isempty(ej)
        shift!(ej)
        # P = P + sparse(ei[1:end-bm],ej,1,nU,nV)
        append!(U1,ei[1:end-bm])
        append!(V1,ej)
      end
    end
  end

  all_matches = [U1 V1]
  unique_matches = unique(all_matches,1)
  U1unique = unique_matches[:,1]
  V1unique = unique_matches[:,2]
  uo = U[U1unique,:]
  vo = V[V1unique,:]
  weights = vec(sum(uo.*vo,2))
  X = sparse(U1unique,V1unique,weights,nU,nV)

  return X
end
