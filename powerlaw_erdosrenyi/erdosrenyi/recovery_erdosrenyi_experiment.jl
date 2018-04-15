include("../../all_code/include_all.jl")
include(joinpath(Pkg.dir(),"handy_julia","handy_julia.jl"))
using NetworkAlign
using MatrixNetworks

function run_erdosrenyi_experiments_keep(iters,nbtrials,bmatch)
n = 50
Pe = vcat(collect(0.0:0.001:0.01),collect(0.02:0.01:0.05))
network_density = [0.1,0.4]
# network_density = [4,6]
noise_model = [1,2]
recov = zeros(length(Pe),length(network_density),length(noise_model),3)

recovEA = zeros(length(Pe),length(network_density),length(noise_model),nbtrials)
recovLR = zeros(length(Pe),length(network_density),length(noise_model),nbtrials)

all_times = zeros(length(Pe),length(network_density),length(noise_model),3)
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
        tic();ma,mb = EigenAlign_while_resid(A,B,s1,s2,s3,1e-12);t = toq()
        r = evaluate_erdosreyni_experiment(A,B,ma,mb)
        all_times[ii,jj,kk,1] += t
        recovEA[ii,jj,kk,trials] = r

        #lowrank_balanced_union_lowrank_evaluation
        #lowrank_balanced_union
        # 2
        if kk == 1
          tic(); ma,mb,D = align_networks_eigenalign(A,B,iters,"lowrank_svd_union",bmatch);t = toq()
        end
        tic(); ma,mb,D = align_networks_eigenalign(A,B,iters,"lowrank_svd_union",bmatch);t = toq()
        r = evaluate_erdosreyni_experiment(A,B,ma,mb)
        all_times[ii,jj,kk,2] += t
        recovLR[ii,jj,kk,trials] = r
        Psprarsity[ii,jj,kk] += D

      end
    end
  end
end
return recovEA,recovLR,Dvals,Psprarsity,all_times
end
