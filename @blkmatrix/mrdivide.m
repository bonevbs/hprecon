function C = mrdivide(A, B)

%C = ( A' \ B' )';

if size(B,1) == size(B,2) && size(A,2) ~= size(B,1)
  error("Incorrect dimensions. Only square matrices are supported.");
end

if isa(A,'blkmatrix')
  if isscalar(B) && all(size(B) == 1)
    C = A * (1 / B);
  elseif isa(B, 'blkmatrix')
    if size(A.A11,2) == size(B.A11,1)
      C = blkmatrix();
      [C.A11, C.A12] = blkmatrix_solve_left(B, A.A11, A.A12);
      [C.A21, C.A22] = blkmatrix_solve_left(B, A.A21, A.A22);
    else
      warning("Blocking does not align. Using full matrix - this might be inefficient.")
      C = full(A)/B;
    end
  else
    C = full(A)/B;
  end
else
  if size(B.A11,1) ~= size(B.A11,2)
    warning("Top left block is not square. Using full matrix.")
    C = A / full(B);
  elseif isscalar(B) && all(size(B) == 1)
    C = A * inv(B);
  else
    [C1, C2] = blkmatrix_solve_left(B, A(:,1:size(B.A11,1)), A(:,size(B.A11,1)+1:end));
    C = [C1, C2];
  end
end


end