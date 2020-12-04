function matching_with_prior(U,V,L)
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
    println("CURRENT COL IS $i")
    ui = U_weights[:,i]
    vi = V_weights[:,i]
    lastid_ui = findfirst(ui.<0)
    lastid_vi = findfirst(vi.<0)

    if lastid_ui == 0 && lastid_vi == 0
      uipositive = U_sortperm[:,i]
      uinegative = []
      vipositive = V_sortperm[:,i]
      vinegative = []

    elseif lastid_vi == 0
      uipositive = U_sortperm[1:lastid_ui-1,i]
      uinegative = []
      vipositive = V_sortperm[:,i]
      vinegative = []

    elseif lastid_ui == 0
      uipositive = U_sortperm[:,i]
      uinegative = []
      vipositive = V_sortperm[1:lastid_vi-1,i]
      vinegative = []

    else
      uipositive = U_sortperm[1:lastid_ui-1,i]
      uinegative = U_sortperm[collect(end:-1:lastid_ui),i]
      vipositive = V_sortperm[1:lastid_vi-1,i]
      vinegative = V_sortperm[collect(end:-1:lastid_vi),i]
    end

    ei1,ej1 = generate_local_matching(uipositive,vipositive,L)
    ei2,ej2 = generate_local_matching(uinegative,vinegative,L)

    ei = vcat(ei1,ei2)
    ej = vcat(ej1,ej2)
    P = P + sparse(ei,ej,1,nU,nV)

  end
  return P
end
function generate_local_matching(ei,ej,L)
  i = 1
  j = 1
  ma = Int64[]
  mb = Int64[]
  nA = length(ei)
  nB = length(ej)

  while i<=nA && j <= nB
    # println("current i is $i and current j is $j")
    a = ei[i]
    b = ej[j]
    # while a==0 && i<nA
    #   i+=1
    #   a=ei[i]
    # end
    # while b==0 && j<nB
    #   j+=1
    #   b=ej[j]
    # end
    #
    if a==0 || b==0
      break;
    end

    if L[a,b] != 0
      append!(ma,a)
      append!(mb,b)
      ei[i] = 0
      ej[j] = 0
    else
      #search for a pair for i
      jinternal = 1
      b = ej[jinternal]
      while jinternal < nB && (b == 0 || L[a,b] == 0)
        jinternal += 1
        b = ej[jinternal]
      end
      #add a,b
      if jinternal <= nB
        if a!=0 && b!=0
          append!(ma,a)
          append!(mb,b)
          ei[i] = 0
          ej[jinternal] = 0
        end
      end

      #search for a pair for j
      iinternal = 1
      a = ei[iinternal]
      b = ej[j]
      while iinternal < nA && (a == 0 || L[a,b] == 0)
        iinternal += 1
        a = ei[iinternal]
      end
      #add a,b
      if iinternal <= nA
        if a!=0 && b!=0
          append!(ma,a)
          append!(mb,b)
          ei[iinternal] = 0
          ej[j] = 0
        end
      end
    end
      #update i
      i += 1
      while i <= nA && ei[i] == 0
        i+=1
      end

      #update j
      j += 1
      while j <= nB && ej[j] == 0
        j += 1
      end
  end
  return ma,mb
end
