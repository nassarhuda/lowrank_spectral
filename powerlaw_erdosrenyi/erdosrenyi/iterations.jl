include("../../all_code/include_all.jl")
include(joinpath(Pkg.dir(),"handy_julia","handy_julia.jl"))
using MatrixNetworks
using Plots
font = Plots.font("sans-serif", 13)
pyplot(size=(400,400),guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)
halfscale = 2
itersvec = 1:15

f = plot(itersvec,ones(length(itersvec)),label="",linewidth=halfscale,ylim=(0,1.2))
g = plot(itersvec,ones(length(itersvec)),label="",linewidth=halfscale,ylim=(0,1.2))


N = [500,1000,1500]
pvals = 20./N
trialnb = 10 #change as needed
rec = zeros(length(N))
mat = zeros(length(N))

recov = zeros(length(itersvec),length(N))
match = zeros(length(itersvec),length(N))

for ni = 1:length(N)
  n = N[ni]
  p = pvals[ni]
  for t = 1:trialnb
    @show (n,t)
    A,B = generate_erdos_reyni_pair(n,p,2,0.5/n)
    s1,s2,s3 = find_parameters(A,B)
    tic(); ei,ej = EigenAlign_while_resid(A,B,s1,s2,s3,1e-12);t = toq()
    rec[ni] += evaluate_erdosreyni_experiment(A,B,ei,ej)
    Ai = A[ei,ei]
    Bi = B[ej,ej]
    mat[ni] += nnz(Ai.*Bi)

    for i = 1:length(itersvec)
      iters = itersvec[i]
      tic(); ma,mb,D = align_networks_eigenalign(A,B,iters,"lowrank_svd_union",3);t = toq()
      r = evaluate_erdosreyni_experiment(A,B,ma,mb)
      @show r
      recov[i,ni] += r
      Ai = A[ma,ma]
      Bi = B[mb,mb]
      match[i,ni] += nnz(Ai.*Bi)
    end
  end
  labval = join(["n = " n])
  plot!(f,itersvec, match[:,ni]./mat[ni],label = labval)
  plot!(g,itersvec, recov[:,ni]./rec[ni],label = labval)
  savefig(f,"itersfig_overlap.pdf")
  savefig(g,"itersfig_recov.pdf")
end