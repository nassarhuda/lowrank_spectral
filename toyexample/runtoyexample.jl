# toy example
include("../all_code/include_all.jl")
A,B = generate_erdos_reyni_pair(50,0.4,1,0.0) # generate two networks to align.
ma,mb = align_networks_eigenalign(A,B,iters,"lowrank_svd_union",3) # alignment step.
# recovery of correctly aligned nodes (as defined in paper):
recovery = evaluate_erdosreyni_experiment(A,B,ma,mb)
