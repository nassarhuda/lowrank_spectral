include("../../../all_code/include_all.jl")
include(joinpath(Pkg.dir(),"handy_julia","handy_julia.jl"))
using NetworkAlign
using MatrixNetworks
# using PyPlot

#get data first

allfiles = ["yeast0_Y2H1.gw",
            "yeast5_Y2H1.gw",
            "yeast10_Y2H1.gw",
            "yeast15_Y2H1.gw",
            "yeast20_Y2H1.gw",
            "yeast25_Y2H1.gw"];

GS3_Netalignbp = zeros(6,6) #
GS3_EigenAlign = zeros(6,6) #
GS3_Netalignmr = zeros(6,6) #

GS3_Netalignbp_lr = zeros(6,6)
GS3_EigenAlign_lr = zeros(6,6) #
GS3_Netalignmr_lr = zeros(6,6) #

time_EA = zeros(6,6) #
time_BP = zeros(6,6) #
time_MR = zeros(6,6) #

time_EA_lr = zeros(6,6) #
time_BP_lr = zeros(6,6)
time_MR_lr = zeros(6,6) #

mharmonic1 = zeros(6,6) #
mharmonic2 = zeros(6,6) #
mharmonic3 = zeros(6,6) #

mharmonic_EA = zeros(6,6) #
mharmonicMR_before = zeros(6,6) #
mharmonicBP_before = zeros(6,6) #

for i = 1:6
  for j = i:6
    println("[i,j] = [$i,$j]")
    A,B,Pref = create_yeast_networks2(allfiles[i],allfiles[j])
    nA = size(A,1)
    nB = size(B,1)

    iters = 8
    ################################################################
    tic()
    L = ones(size(A,1),size(B,1))
    L = sparse(L)
    S,w,li,lj = netalign_setup(A,B,L)
    # align networks
    a = 1;
    b = 1;
    xbest1,st,status,hist = netalignmr(S,w,1,1,li,lj,0.4,25,1,50,true)
    ma,mb = edge_list(bipartite_matching(xbest1,li,lj));
    time_MR[i,j] = toq()
    GS3_Netalignmr[i,j] = gs3_evaluate(A,B,ma,mb)
    Presult_bp = sparse(ma,mb,1,nA,nB)
    mharmonicMR_before[i,j] = mean_harmonic(Presult_bp,Pref)
    ################################################################
    tic()
    L = ones(size(A,1),size(B,1))
    L = sparse(L)
    S,w,li,lj = netalign_setup(A,B,L)
    # align networks
    a = 1;
    b = 1;
    stepm = 25;
    rtype = 1;
    maxiter = 5;
    verbose = true;
    gamma = 0.4;
    xbest3,flag,reshist = netalignbp(S,w,a,b,li,lj,0.99,1,maxiter,true)
    ma,mb = edge_list(bipartite_matching(xbest3,li,lj));
    time_BP[i,j] = toq()
    GS3_Netalignbp[i,j] = gs3_evaluate(A,B,ma,mb)
    Presult_bp = sparse(ma,mb,1,nA,nB)
    mharmonicBP_before[i,j] = mean_harmonic(Presult_bp,Pref)
    ################################################################
    tic(); ma,mb,D = align_networks_eigenalign(A,B,iters,"lowrank_svd_union",3); t = toc()
    GS3_EigenAlign_lr[i,j] = gs3_evaluate(A,B,ma,mb)
    time_EA_lr[i,j] = t
    Presult_bp = sparse(ma,mb,1,nA,nB)
    mharmonic1[i,j] = mean_harmonic(Presult_bp,Pref)


    # ################################################################
    tic()
    s1,s2,s3 = find_parameters(A,B)
    c1 = s1+s2-2s3; c2 = s3-s2; c3 = s2
    Uk,Vk,Wk,W1,W2 = decomposeX_balance_allfactors(A,B,8,c1,c2,c3)
    Un,Vn = split_balanced_decomposition(Uk,Wk,Vk)
    Xout = newbound_rounding_lowrank_evaluation_relaxed(Un,Vn,3)
    S,w,li,lj = netalign_setup(A,B,Xout)
    # align networks
    a = 1;
    b = 1;
    stepm = 25;
    rtype = 1;
    maxiter = 50;
    verbose = true;
    gamma = 0.4;
    xbest1,st,status,hist = netalignmr(S,w,a,b,li,lj,gamma,stepm,rtype,maxiter,verbose)
    ma,mb = edge_list(bipartite_matching(xbest1,li,lj));
    time_MR_lr[i,j] = toq()
    GS3_Netalignmr_lr[i,j] = gs3_evaluate(A,B,ma,mb)
    Presult_bp = sparse(ma,mb,1,nA,nB)
    mharmonic2[i,j] = mean_harmonic(Presult_bp,Pref)

    # ################################################################
    tic()
    s1,s2,s3 = find_parameters(A,B)
    c1 = s1+s2-2s3; c2 = s3-s2; c3 = s2
    Uk,Vk,Wk,W1,W2 = decomposeX_balance_allfactors(A,B,8,c1,c2,c3)
    Un,Vn = split_balanced_decomposition(Uk,Wk,Vk)
    Xout = newbound_rounding_lowrank_evaluation_relaxed(Un,Vn,3)
    S,w,li,lj = netalign_setup(A,B,Xout)
    # align networks
    a = 1;
    b = 1;
    stepm = 25;
    rtype = 1;
    maxiter = 50;
    verbose = true;
    gamma = 0.4;
    xbest3,flag,reshist = netalignbp(S,w,a,b,li,lj,0.99,1,maxiter,true)
    ma,mb = edge_list(bipartite_matching(xbest3,li,lj));
    time_BP_lr[i,j] = toq()
    GS3_Netalignbp_lr[i,j] = gs3_evaluate(A,B,ma,mb)
    Presult_bp = sparse(ma,mb,1,nA,nB)
    mharmonic3[i,j] = mean_harmonic(Presult_bp,Pref)

    #################################################################

    tic(); ma,mb = EigenAlign_while_resid(A,B,s1,s2,s3,1e-12); t= toc()
    GS3_EigenAlign[i,j] = gs3_evaluate(A,B,ma,mb)
    time_EA[i,j] = t
    Presult_bp = sparse(ma,mb,1,nA,nB)
    mharmonic_EA[i,j] = mean_harmonic(Presult_bp,Pref)

  end
end
