include("normout.jl")

#TODO: define input types
#TODO:

## julia code for NSD:
function NSD(A,B,z,w,alpha,iters,m,n)
  #compute A_tilde and B_tilde, no need to transpose (just use A'x later on)
  #A = A./sum(A,2)
  #B = B./sum(B,2)

  s = size(w,2)
  Z = zeros(n,iters)#A
  W = zeros(m,iters)#B
  Sim = zeros(m,n)

  # for i = 1:s
  # initially this for loop was for multiple low rank factors
  # i removed it because I just wanted quick access to the lowrank factors
  i = 1
    @show s
    W[:,1] = w
    Z[:,1] = z
    temp = w
    f = alpha
    for k = 2:iters
      temp = B'*temp
      W[:,k] = f*temp
      Z[:,k] = A'*Z[:,k-1]
      f = f*alpha
    end
    w_last = f*B'*temp
    z_last = A'*Z[:,iters]

    X = (1-alpha)*W*Z' +w_last*z_last'

    Sim = Sim + X
  # end
  # added by me on a second phase
  Wtoreturn = hcat((1-alpha)*W,w_last)
  Ztoreturn = hcat(Z,z_last)
  return Sim,Wtoreturn,Ztoreturn
end


##############

function NSD_setup(A, B, preiters, iters, alpha)
# preiters sets up the inputs wi and zi
# iters is n in the paper

# keep same naming convention they used!

n = size(A,1)
vecA = ones(Float64,n) ./ n

m = size(B,1)
vecB = ones(Float64,m) ./ m

# set up A and B
# @time A = sparse(max(A,A'))
# @time B = sparse(max(B,B'))
@time A += A'
@time B += B'
P = normout(A)
Q = normout(B)

vecsA = zeros(n,preiters+1)
vecsA[:,1] = vecA

vecsB = zeros(m,preiters+1)
vecsB[:,1] = vecB

# compute vecsA and vecsB
for i = 1: preiters
    vecA = P'*vecA
    vecB = Q'*vecB
    vecsA[:,i+1] = vecA
    vecsB[:,i+1] = vecB
end

ivecsA = vecsA[:,preiters+1]
ivecsB = vecsB[:,preiters+1]

@time X,Wtoreturn,Ztoreturn = NSD(P,Q,ivecsA,ivecsB,alpha,iters,m,n)
return X,Wtoreturn,Ztoreturn
end
