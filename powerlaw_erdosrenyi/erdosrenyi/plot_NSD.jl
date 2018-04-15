nbtrials = 100
iters = 8
bmatch = 3
# include("run_erdosrenyi_experiments_withNSD.jl")
# recov,Dvals,Psprarsity,all_times = run_erdosrenyi_experiments_withNSD(iters,nbtrials,bmatch)

using Plots
font = Plots.font("sans-serif")
pyplot(size=(300,300),guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)
halfscale = 1.5

Pe = vcat(collect(0.0:0.001:0.01),collect(0.02:0.01:0.05))
network_density = [4,6]
noise_model = [1,2]

recovreal = recov./nbtrials
labels = ["NSD+GM","NSD+LR","EA+LR","EA+GM"]
# plot:
plot()
for i = 1:length(network_density)
  for j = 1:length(noise_model)
    for k = 1:4
      ea = vec(recovreal[:,i,j,k])
      # plot!(Pe,ea,label=labels[k])
      plot!(Pe,ea,label=labels[k], linewidth=halfscale)
    end
    xlabel!("noise level")
    # ylabel!("recovery rate of true mappings")
    # title!("density = $(network_density[i]), noise model = $(noise_model[j])")
    # imgfile = join(["_bval_",bmatch,"_eff_newbound_recovery_analysis", network_density[i], "_noise", noise_model[j],"_",iters,"iters.pdf"])
    imgfile = join(["EA_NSD_",i,"_",j,".pdf"])
    savefig(imgfile)
    close()
    plot()
  end
end
