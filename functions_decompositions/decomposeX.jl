function decomposeX_balance_allfactors(A,B,k,c1,c2,c3)
  nA = size(A,1)
  nB = size(B,1)

  u = ones(nA,1)./nA
  v = ones(nB,1)./nB
  du = zeros(k)
  dv = zeros(k)
  du[end] = 1
  dv[end] = 1

  eA = ones(nA,1)
  eB = ones(nB,1)

  U = zeros(nA,k)
  U[:,end] = eA
  for i = k-1:-1:2
    U[:,i] = A*U[:,i+1]
    mui = maximum(U[:,i])
    U[:,i] = U[:,i]./mui
    du[i] = du[i+1]*mui
  end
  U[:,1] = u
  rksums = zeros(k-1)
  du[1] = 1
  for i = 1:k-1
    rksums[i] = du[1]*sum(U[:,1])
    U[:,1] = A*U[:,1]
    mui = maximum(U[:,1])
    U[:,1] = U[:,1]./mui
    du[1] = du[1]*mui
  end
  usums = vec(sum(U,1))
  usums = usums.*du
  # println("done with U")

  V = zeros(nB,k)
  V[:,end] = eB
  for i = k-1:-1:2
    V[:,i] = B*V[:,i+1]
    mvi = maximum(V[:,i])
    V[:,i] = V[:,i]./mvi
    dv[i] = dv[i+1]*mvi
  end
  V[:,1] = v
  hksums = zeros(k-1)
  dv[1] = 1
  for i = 1:k-1
    hksums[i] = dv[1]*sum(V[:,1])
    V[:,1] = B*V[:,1]
    mvi = maximum(V[:,1])
    V[:,1] = V[:,1]./mvi
    dv[1] = dv[1]*mvi
  end
  vsums = vec(sum(V,1))
  vsums = vsums.*dv

  rhos = copy(du)
  gams = copy(dv)

    # create W
    W_prev = 1*rhos[1]*gams[1]
    W1 = 1
    W2 = 1

    scalerow = 1
    scalecol = 1

    for i = 2:k
      localk = i-1
      rk = zeros(localk)
      rk[1] = rksums[localk]
      lt = length(2:localk)
      rk[2:localk] = usums[end-lt+1:end]

      hk = zeros(localk)
      hk[1] = hksums[localk]
      hk[2:localk] = vsums[end-lt+1:end]

      Du = Diagonal(rhos[i]./rhos[1:i-1])
      Dv = Diagonal(gams[i]./gams[1:i-1])

      v1 = Dv*hk
      v2 = Du*rk

      W_curr = vcat(hcat(c1*W_prev,scalecol*c2*W_prev*v1),
               hcat(scalerow*c2*v2'*W_prev,scalerow*scalecol*c3*v2'*W_prev*v1))

#      W_curr = vcat(hcat(c1*W_prev,scalecol*c2*W_prev*Dv*hk),
#                    hcat(scalerow*c2*rk'*Du*W_prev,scalerow*scalecol*c3*rk'*Du*W_prev*Dv*hk))

      W1 = vcat(hcat(W_prev,zeros(size(W_prev,1),localk)),
                  hcat(zeros(1,size(W_prev,2)),rk'))
      W2 = vcat(hcat(c1*eye(size(W_prev,1)),c2*hk),
                  hcat(c2*W_prev,c3*W_prev*hk))

      W_curr = W_curr/W_curr[end,end]

      W_prev = W_curr
    end
  return U,V,W_prev,W1,W2
end

# old version of decomposeX
function decomposeX(A,B,k,c1,c2,c3)
  println("decomposeX")
  nA = size(A,1)
  nB = size(B,1)

  u = ones(nA,1)./nA
  v = ones(nB,1)./nB

  eA = ones(nA,1)
  eB = ones(nB,1)

  U = zeros(nA,k)
  U[:,end] = eA
  for i = k-1:-1:2
    U[:,i] = A*U[:,i+1]
  end
  U[:,1] = u
  rksums = zeros(k-1)
  for i = 1:k-1
    rksums[i] = sum(U[:,1])
    U[:,1] = A*U[:,1]
  end
  usums = sum(U,1)
# println("done with U")
  V = zeros(nB,k)
  V[:,end] = eB
  for i = k-1:-1:2
    V[:,i] = B*V[:,i+1]
  end
  V[:,1] = v
  hksums = zeros(k)
  for i = 1:k-1
    hksums[i] = sum(V[:,1])
    V[:,1] = B*V[:,1]
  end
  vsums = sum(V,1)
# println("done with U")
  # create W
  W_prev = 1
  W1 = 1
  W2 = 1
  @show usums
  @show vsums

  for i = 2:k
# println("i is $i out of $k")
    localk = i-1
    rk = zeros(localk)
    rk[1] = rksums[localk]
    lt = length(2:localk)
    rk[2:localk] = usums[end-lt+1:end]
    @show rk

    # rk[1] = sum(A^(localk-1)*u)
    # for j = 2:localk
    #   rk[j] = sum(A^(localk-j)*eA)
    # end


    hk = zeros(localk)
    hk[1] = hksums[localk]
    hk[2:localk] = vsums[end-lt+1:end]
    @show hk
    # hk[1] = sum(B^(localk-1)*v)
    # for j = 2:localk
    #   hk[j] = sum(B^(localk-j)*eB)
    # end

    # println("Wcalc start")
    W_curr = vcat(hcat(c1*W_prev,c2*W_prev*hk),
                  hcat(c2*rk'*W_prev,c3*rk'*W_prev*hk))

    W1 = vcat(hcat(W_prev,zeros(size(W_prev,1),localk)),
                hcat(zeros(1,size(W_prev,2)),rk'))
    W2 = vcat(hcat(c1*eye(size(W_prev,1)),c2*hk),
                hcat(c2*W_prev,c3*W_prev*hk))

    W_prev = W_curr
    # println("Wcalc end")
  end
  return U,V,W_prev,W1,W2
end

function decomposeX_normalized(A,B,k,c1,c2,c3)
  nA = size(A,1)
  nB = size(B,1)

  u = ones(nA,1)./nA
  v = ones(nB,1)./nB

  eA = ones(nA,1)
  eB = ones(nB,1)

  U = zeros(nA,k)
  U[:,end] = eA
  U[:,end] = U[:,end]./sum(U[:,end])
  for i = k-1:-1:2
    U[:,i] = A*U[:,i+1]
    U[:,i] = U[:,i]./sum(U[:,i])
  end
  U[:,1] = u
  rksums = zeros(k-1)
  for i = 1:k-1
    rksums[i] = sum(U[:,1])
    U[:,1] = A*U[:,1]
  end
  U[:,1] = U[:,1]./sum(U[:,1])
  usums = sum(U,1)
# println("done with U")
  V = zeros(nB,k)
  V[:,end] = eB
  V[:,end] = V[:,end]./sum(V[:,end])
  for i = k-1:-1:2
    V[:,i] = B*V[:,i+1]
    V[:,i] = V[:,i]./sum(V[:,i])
  end
  V[:,1] = v
  hksums = zeros(k)
  for i = 1:k-1
    hksums[i] = sum(V[:,1])
    V[:,1] = B*V[:,1]
  end
  V[:,1] = V[:,1]./sum(V[:,1])
  vsums = sum(V,1)
# println("done with U")
  # create W
  W_prev = 1
  W1 = 1
  W2 = 1
  @show usums
  @show vsums
  @show hksums
  @show rksums
  for i = 2:k
# println("i is $i out of $k")
    localk = i-1
    rk = zeros(localk)
    rk[1] = rksums[localk]
    lt = length(2:localk)
    rk[2:localk] = usums[end-lt+1:end]
    # rk[1] = sum(A^(localk-1)*u)
    # for j = 2:localk
    #   rk[j] = sum(A^(localk-j)*eA)
    # end


    hk = zeros(localk)
    hk[1] = hksums[localk]
    hk[2:localk] = vsums[end-lt+1:end]
    # hk[1] = sum(B^(localk-1)*v)
    # for j = 2:localk
    #   hk[j] = sum(B^(localk-j)*eB)
    # end

    # println("Wcalc start")
    W_curr = vcat(hcat(c1*W_prev,c2*W_prev*hk),
                  hcat(c2*rk'*W_prev,c3*rk'*W_prev*hk))

    W1 = vcat(hcat(W_prev,zeros(size(W_prev,1),localk)),
                hcat(zeros(1,size(W_prev,2)),rk'))
    W2 = vcat(hcat(c1*eye(size(W_prev,1)),c2*hk),
                hcat(c2*W_prev,c3*W_prev*hk))

    W_prev = W_curr
    # println("Wcalc end")
  end
  return U,V,W_prev,W1,W2
end

function decomposeX_with_uv(A,B,k,c1,c2,c3,u,v)
  nA = size(A,1)
  nB = size(B,1)

  # u = ones(nA,1)./nA
  # v = ones(nB,1)./nB

  eA = ones(nA,1)
  eB = ones(nB,1)

  U = zeros(nA,k)
  U[:,end] = eA
  for i = k-1:-1:2
    U[:,i] = A*U[:,i+1]
  end
  U[:,1] = u
  rksums = zeros(k-1)
  for i = 1:k-1
    rksums[i] = sum(U[:,1])
    U[:,1] = A*U[:,1]
  end
  usums = sum(U,1)
# println("done with U")
  V = zeros(nB,k)
  V[:,end] = eB
  for i = k-1:-1:2
    V[:,i] = B*V[:,i+1]
  end
  V[:,1] = v
  hksums = zeros(k)
  for i = 1:k-1
    hksums[i] = sum(V[:,1])
    V[:,1] = B*V[:,1]
  end
  vsums = sum(V,1)
# println("done with U")
  # create W
  W_prev = 1
  W1 = 1
  W2 = 1
  for i = 2:k
# println("i is $i out of $k")
    localk = i-1
    rk = zeros(localk)
    rk[1] = rksums[localk]
    lt = length(2:localk)
    rk[2:localk] = usums[end-lt+1:end]
    # rk[1] = sum(A^(localk-1)*u)
    # for j = 2:localk
    #   rk[j] = sum(A^(localk-j)*eA)
    # end


    hk = zeros(localk)
    hk[1] = hksums[localk]
    hk[2:localk] = vsums[end-lt+1:end]
    # hk[1] = sum(B^(localk-1)*v)
    # for j = 2:localk
    #   hk[j] = sum(B^(localk-j)*eB)
    # end

    # println("Wcalc start")
    W_curr = vcat(hcat(c1*W_prev,c2*W_prev*hk),
                  hcat(c2*rk'*W_prev,c3*rk'*W_prev*hk))

    W1 = vcat(hcat(W_prev,zeros(size(W_prev,1),localk)),
                hcat(zeros(1,size(W_prev,2)),rk'))
    W2 = vcat(hcat(c1*eye(size(W_prev,1)),c2*hk),
                hcat(c2*W_prev,c3*W_prev*hk))

    W_prev = W_curr
    # println("Wcalc end")
  end
  return U,V,W_prev,W1,W2
end
