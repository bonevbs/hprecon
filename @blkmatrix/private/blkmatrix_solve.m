function [x1, x2] = blkmatrix_solve(A, b1, b2)
  x1 = A.A11\b1;
  x2 = b2 - A.A21*x1;
  if isempty(A.S22)
    if isa(A.A11, 'hss') && isa(A.A22, 'hss')
      % this should be replaced by a method that uses randomized sampling
      A.S22 = hss(A.A22 - A.A21*(A.A11\A.A12));
    else
      A.S22 = A.A22 - A.A21*(A.A11\A.A12);
    end
  end
  x2 = A.S22\x2;
  x1 = x1 - A.A11\(A.A12*x2);
end