
function mean_harmonic(P1,P2)
  x1 = nnz(P1.*P2)/nnz(P2)
  x2 = nnz(P1.*P2)/nnz(P1)
  mharmonic = 2/((1/x1)+(1/x2))
  return mharmonic
end
