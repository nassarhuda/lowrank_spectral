using Plots

function find_percentile(v, pct)
           len = length(v)
           ind = floor(Int64,pct/100*len)
           newarr = sort(v);
           val = newarr[ind];
           return val
end

recovEA,recovLR,Dvals,Psprarsity,all_times = run_powerlaw_experiments_keep(8,200,3);
font = Plots.font("sans-serif", 13)
pyplot(size=(400,400),guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)
halfscale = 1
Pe = vcat(collect(0.0:0.001:0.01),collect(0.02:0.01:0.05))
network_density = [4,6]
noise_model = [1,2]
r20 = zeros(15)
r50 = zeros(15)
r80 = zeros(15)
for jj = 1:length(network_density)
  for kk = 1:2 # noise model 2 only
    r1 = recovEA[:,jj,kk,1:100]
    for ii = 1:length(Pe)
      r20[ii] = find_percentile(r1[ii,:],20)
      r50[ii] = find_percentile(r1[ii,:],50)
      r80[ii] = find_percentile(r1[ii,:],80)
    end
    plot(Pe,r50,label="EigenAlign", linewidth=2*halfscale, color = 1,yticks=0.5:0.25:1)
    plot!(Pe,r20, label="", linestyle = :dash, linewidth=halfscale, color = 1)
    plot!(Pe,r80, label="", linestyle = :dash, linewidth=halfscale, color = 1)

    r1 = recovLR[:,jj,kk,1:100]
    for ii = 1:length(Pe)
      r50[ii] = find_percentile(r1[ii,:],50)
      r20[ii] = find_percentile(r1[ii,:],20)
      r80[ii] = find_percentile(r1[ii,:],80)
    end
    plot!(Pe,r50,label="low rank ", linewidth=2*halfscale, color = 2)
    plot!(Pe,r20, label="", linestyle = :dash, linewidth=halfscale, color = 2)
    plot!(Pe,r80, label="", linestyle = :dash, linewidth=halfscale, color = 2)
    xlabel!("noise level")
    imgfile = join(["powerlaw",network_density[jj],".pdf"])
    savefig(imgfile)
    close()
  end
end
