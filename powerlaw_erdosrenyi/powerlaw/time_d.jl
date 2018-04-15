include("../../all_code/include_all.jl")
include(joinpath(Pkg.dir(),"handy_julia","handy_julia.jl"))

trialnumbers = 10

N = [50, 100, 200, 500, 1000, 2000, 5000, 10000]

time1 = zeros(length(N),trialnumbers)
time2 = zeros(length(N))
time3 = zeros(length(N),trialnumbers)
timen = zeros(length(N))
Dvals = zeros(length(N))
Psparsity = zeros(length(N))

matchingtime_EA = zeros(length(N),trialnumbers)
matchingtime_LR = zeros(length(N),trialnumbers)
#

recov1_store = zeros(length(N),trialnumbers)
recov3_store = zeros(length(N),trialnumbers)
DVALS = zeros(length(N),trialnumbers)
DALL = zeros(length(N)*trialnumbers,3)
Psparsityn = zeros(length(N),trialnumbers)

# for precompiling
A,B = generate_powerlaw_pair(50,4,2,0.1)
s1,s2,s3 = find_parameters(A,B)

iters = 8
bm = 3
k = 0
for i = 1:length(N)
  @show i
  n = N[i]
  @show n

  Pe = 0.5/n

  if i == 1 #precompile
    ma,mb,D,tt = align_networks_eigenalign(A,B,iters,"lowrank_svd_union",bm);
    ei,ej = EigenAlign_while_resid(A,B,s1,s2,s3,1e-12)
  end
  for trialnb = 1:trialnumbers
  k += 1

    A,B = generate_powerlaw_pair(n,4,2,Pe)
    s1,s2,s3 = find_parameters(A,B)

    tic(); ma,mb,D,timematching = align_networks_eigenalign(A,B,iters,"lowrank_svd_union",bm);t = toq()
    matchingtime_LR[i,trialnb] = timematching
    time1[i,trialnb] = t

    if n <= 1000
      tic();ei,ej,matchingtime = EigenAlign_while_resid(A,B,s1,s2,s3,1e-12);t = toq()
      matchingtime_EA[i,trialnb] += matchingtime
      time3[i,trialnb] = t
      
      #all things:
      tic(); ma,mb,D,timematching = align_networks_eigenalign(A,B,iters,"lowrank_balanced_best_eff",0);t = toc()

      DALL[k,:] = D
    end

  end
end



#### plotting:
ids = 1:5
time1avg = mean(time1,2)
matchingtime_LR_avg = mean(matchingtime_LR,2)
time3avg = mean(time3,2)
matchingtime_EA_avg = mean(matchingtime_EA,2)
NN = N
plot(NN,time1avg./ntrials,label="low rank", linewidth=2*halfscale,color = 2)
plot!(NN,matchingtime_LR_avg./ntrials,label="", linewidth=halfscale, color = 2,linestyle = :dash)
plot!(yticks=(vec(1:length(1e2))),NN[ids],time3avg[ids]./ntrials,label="EigenAlign", linewidth=2*halfscale, color = 1)
plot!(NN[ids],matchingtime_EA_avg[ids]./ntrials,label="", linewidth=halfscale,color = 1,linestyle = :dash)
plot!(yticks=[0.001,0.01,0.1,1,10,100])
yaxis!("Time (sec)",:log10)
xaxis!("size of the network",:log10)
imgfile = join(["218_218_time.pdf"])
savefig(imgfile)
close()

dvals = DALL[1:50,:]
NN = [50, 100, 200, 500, 1000]

f = plot()
for i = 1:5
  dv = dvals[(i-1)*10+1 : i*10,1]
  scatter!(f,ones(Int,10)*NN[i],dv,color= 2,markersize= 6)#,label="",markersize= 10,color=2)#markershape=:star)
end
xaxis!(f,:log10)
title!(f,"computed approximation bound D")

g = plot()
for i = 1:5
  dv = dvals[(i-1)*10+1 : i*10,2]
  scatter!(g,ones(Int,10)*N[i],dv,color = 1,markersize= 6)#,label="",markersize= 10,color=2)#markershape=:star)
end
title!(g,"actual LR approximation value")
xaxis!(g,:log10)

h = plot()
for i = 1:5
  dv = dvals[(i-1)*10+1 : i*10,3]
  scatter!(h,ones(Int,10)*N[i],dv,color = 3,markersize= 6)#,label="",markersize= 10,color=2)#markershape=:star)
end
title!(h,"actual GM approximation value")
xaxis!(h,:log10)

l = @layout([a;b;c])
p = plot(
      f,
      g,
      h,
      leg=false,
      layout=l)

savefig("218_218_dvals.pdf")
close()
