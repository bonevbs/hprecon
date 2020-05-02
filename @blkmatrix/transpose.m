function B = transpose(A)
  B = blkmatrix(A.A11.', A.A21.', A.A12.', A.A22.');
end