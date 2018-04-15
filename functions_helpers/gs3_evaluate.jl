function gs3_evaluate(A,B,mi,mj)
  mi = collect(mi)
  mj = collect(mj)
  num_matches = length(mi)
  NCV = (length(unique(mi)) + length(unique(mj))) / (size(A, 1) + size(B, 1));

  #ii,jj,vv = findnz(A)
  #ei = mi[ii]
  #ej = mi[jj]
  #A_induced = sparse(ei,ej,vv)

  #ii,jj,vv = findnz(B)
  #ei = mj[ii]
  #ej = mj[jj]
  #B_induced = sparse(ei,ej,vv)

  A_induced = A[mi, mi];
  B_induced = B[mj, mj];
  D = A_induced + B_induced;
  D = D - spdiagm(diag(D))
  ii,jj = ind2sub(D,find(D.==2))
  align_graph = sparse(ii, jj, 1, num_matches, num_matches)

  M_align = nnz(align_graph) / 2
  # T_align = full(trace(align_graph^3)/6);

  GS3 = M_align / (nnz(A_induced)/2 + nnz(B_induced)/2 - M_align)
  NCV_GS3 = 1/(((1/GS3)+(1/NCV))/2)
  return NCV_GS3
end
