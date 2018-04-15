include("../../all_code/include_all.jl")
include(joinpath(Pkg.dir(),"handy_julia","handy_julia.jl"))
using NetworkAlign
using MatrixNetworks

function run_erdosrenyi_experiments_withNSD(iters,nbtrials,bmatch)
n = 50
Pe = vcat(collect(0.0:0.001:0.01),collect(0.02:0.01:0.05))
network_density = [0.1,0.4]
# network_density = [4,6]
noise_model = [1,2]
recov = zeros(length(Pe),length(network_density),length(noise_model),4)
all_times = zeros(length(Pe),length(network_density),length(noise_model),4)
Psprarsity = zeros(length(Pe),length(network_density),length(noise_model))
Dvals = zeros(length(Pe),length(network_density),length(noise_model))
for trials = 1:nbtrials
  for ii = 1:length(Pe)
    for jj = 1:length(network_density)
      for kk = 1:length(noise_model)
        println("[i,j,k] = [$ii,$jj,$kk] - trial is $trials")
        A,B = generate_erdos_reyni_pair(n,network_density[jj],noise_model[kk],Pe[ii])
        # A,B = generate_powerlaw_pair(n,network_density[jj],noise_model[kk],Pe[ii])
        s1,s2,s3 = find_parameters(A,B)

        # 1
        # tic();ma,mb = EigenAlign_while_resid(A,B,s1,s2,s3,1e-12);t = toq()
        # tic();ma,mb,Xmat,weight,resid = EigenAlign_nokron(A,B,s1,s2,s3,200);t = toq()

        tic();Xout,Un,Vn = NSD_setup(A, B, 1, iters-1, 0.2)
        ma,mb = greedy_lowrank(Un,Vn);t= toq()

        # tic(); ma,mb = EigenAlign_matching(A,B,s1,s2,s3,iters);t = toq()
        @show size(A)
        @show size(B)
        @show maximum(ma)
        @show maximum(mb)
        @show minimum(mb)
        @show minimum(mb)

        r = evaluate_erdosreyni_experiment(A,B,ma,mb)
        all_times[ii,jj,kk,1] += t
        recov[ii,jj,kk,1] += r

        tic();X = newbound_rounding_lowrank_evaluation_relaxed(Un,Vn,3) #bmatch
        ma,mb = edge_list(bipartite_matching(X));t=toq()
        r = evaluate_erdosreyni_experiment(A,B,ma,mb)
        all_times[ii,jj,kk,2] += t
        recov[ii,jj,kk,2] += r


        #lowrank_balanced_union_lowrank_evaluation
        #lowrank_balanced_union
        # 2
        if kk == 1
          tic(); ma,mb,D = align_networks_eigenalign(A,B,iters,"lowrank_svd_union",bmatch);t = toq()
        end
        tic(); ma,mb,D = align_networks_eigenalign(A,B,iters,"lowrank_svd_union",bmatch);t = toq()
        r = evaluate_erdosreyni_experiment(A,B,ma,mb)
        all_times[ii,jj,kk,3] += t
        recov[ii,jj,kk,3] += r
        Psprarsity[ii,jj,kk] += D

        # 3
        if kk == 1 #force precompile
          tic(); ma,mb,D = align_networks_eigenalign(A,B,iters,"lowrank_balanced_best_eff",bmatch);t = toq()
        end
        tic(); ma,mb,D = align_networks_eigenalign(A,B,iters,"greedy",bmatch);t = toq()
        r = evaluate_erdosreyni_experiment(A,B,ma,mb)
        all_times[ii,jj,kk,4] += t
        recov[ii,jj,kk,4] += r
        Dvals[ii,jj,kk] += D
      end
    end
  end
end
return recov,Dvals,Psprarsity,all_times
end
