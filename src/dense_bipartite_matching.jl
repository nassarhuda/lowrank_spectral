

function bipartite_matching_primal_dual(X::Matrix{Float64};tol::Float64=1e-8,
    normalize_weights::Bool=false)
#to get the access pattern right, we must match the right hand side to the left hand side. 

m,n = size(X)
@assert m >= n  #error occurs when m < n 

# variables used for the primal-dual algorithm
# normalize ai values # updated on 2-19-2019
if normalize_weights
X ./= maximum(abs.(X))
end

#initialize variables
alpha=zeros(Float64,n)
bt=zeros(Float64,m+n)#beta
match1 = zeros(Int64,n)
match2 = zeros(Int64,n+m)
queue=zeros(Int64,n)
t=zeros(Int64,m+n)
tmod = zeros(Int64,m+n)
ntmod=0

# initialize the primal and dual variables

for j = 1:n
    for i=1:m
        if X[i,j] > alpha[j]
            alpha[j]=X[i,j]
        end
    end
end


# dual variables (bt) are initialized to 0 already
# match1 and match2 are both 0, which indicates no matches

j=1
@inbounds while j<=n

    for i=1:ntmod
        t[tmod[i]]=0
    end
    ntmod=0
    # add i to the stack
    head=1
    tail=1
    queue[head]=j
    while head <= tail && match1[j]==0 #queue empty + i is unmatched

        k=queue[head]
        #println("begining of queue loop")

        for i=1:m+1 #iterate over column k

            if i == m+1 #check the dummy node
                i = k+m
            end
            if i == k+m #dummy nodes don't have weight
                if 0.0 < alpha[k] + bt[i] - tol
                    continue
                end
            elseif X[i,k] < alpha[k] + bt[i] - tol
                continue
            end # skip if tight

            if t[i]==0
                tail=tail+1 #put the potential match in the queue
                if tail <= m
                    queue[tail]=match2[i]
                end
                t[i]=k  #try vertex k for vertex j
                ntmod=ntmod+1
                tmod[ntmod]=i
                if match2[i]<1  #if i is unmatched
                    while i>0   #unfurl out to get an augmented path
                        match2[i]=t[i]
                        k=t[i]
                        temp=match1[k]
                        match1[k]=i
                        i=temp
                    end
                break
            end
        end
    end
    head=head+1
    end

    #if node j was unable to be matched, update flows and search for new augmenting path
    if match1[j] < 1
        theta=Inf
        for i=1:head-1
            t1=queue[i]
            for t2=1:m
            if t[t2] == 0 && alpha[t1] + bt[t2] - X[t2,t1] < theta
            theta = alpha[t1] + bt[t2] - X[t2,t1]
            end
            end
            #check t1's dummy node
            if t[t1 + m] == 0 && alpha[t1] + bt[t1 + m] < theta
            theta = alpha[t1] + bt[t1 + m]
            end
            end

            for i=1:head-1
            alpha[queue[i]] -= theta
            end
            for i=1:ntmod
            bt[tmod[i]] += theta
            end
            continue
        end

        j=j+1
    end

    #count
    val=0.0
    for j=1:n
        for i=1:m
            if i==match1[j]
                val=val+X[i,j]
            end
        end
    end

    #count how many are properly matched
    noute = 0
    for j=1:n
        if match1[j]<=m
            noute=noute+1
        end
    end

    return val,noute,match1, match2
end