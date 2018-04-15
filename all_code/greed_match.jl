function greedy_match(X)
# % Input arguments:
# % - X: the matrix with the similarity scores (similarity matrix).
# %     Note that element X(i,j) is the similarity score of node i in B
# %     and node j in A; if B has m nodes and A has n nodes then X is an
# %     m x n matrix.
# % Output arguments:
# % - M: the sparse matrix with the matches: M(i,j) = 1.0 iff node i in B
# %     matches with node j in A. m x n matrix (same dimensions with X).
# % - dt: the time in seconds for the operation.
#
#
# % Giorgos Kollias and Shahin Mohammadi
# % Department of Computer Science, Purdue University

m, n = size(X)
N = m * n
x = X[:]
minSize = min(m, n)
usedRows = zeros(m)
usedCols = zeros(n)

maxList = zeros(minSize)
row = zeros(Int,minSize)
col = zeros(Int,minSize)
# y = sort(x,rev=true)
ix = sortperm(x,rev=true)
y = x[ix]

matched = 1
index = 1
while matched <= minSize
    ipos = ix[index] # position in the original vectorized matrix
    jc = ceil(Int,ipos / m)
    ic = ipos - (jc - 1) * m
    if ic == 0
      ic = 1
    end

    if usedRows[ic] != 1 && usedCols[jc] != 1
      # matched;
      row[matched] = ic
      col[matched] = jc
      maxList[matched] = x[index]
      usedRows[ic] = 1
      usedCols[jc] = 1
      matched += 1
    end
    index += 1
end

data = ones(minSize)
M = sparse(row, col, data, m, n);
return row,col,M
end
