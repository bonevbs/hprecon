function [x1, x2] = blksolve_left(blkA, b1, b2)
  A11 = blkA.A11; A12 = blkA.A12; A21 = blkA.A21; A22 = blkA.A22;
  x1 = b1/A11;
  x2 = b2 - x1*A12;
  S22d = hss(A22 - A21*(A11\A12));
  x2 = x2/S22d;
  x1 = x1 - (x2*A21)/A11;
end