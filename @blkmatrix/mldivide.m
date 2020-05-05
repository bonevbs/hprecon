function C = mldivide(A, B)

if size(A,1) == size(A,2) && size(A,2) ~= size(B,1)
  error("Incorrect dimensions. Only square matrices are supported.");
end

if isa(B,'blkmatrix')
  if isscalar(A) && all(size(A) == 1)
    C = (1 / A) * B;
  elseif isa(A, 'blkmatrix')
    if size(A.A11,2) == size(B.A11,1)
      C = blkmatrix();
      [C.A11, C.A21] = blkmatrix_solve(A, B.A11, B.A21);
      [C.A12, C.A22] = blkmatrix_solve(A, B.A12, B.A22);
    else
      warning("Blocking does not align. Using full matrix - this might be inefficient.")
      C = A\full(B);
    end
  else
    C = A\full(B);
  end
else
  if size(A.A11,1) ~= size(A.A11,2)
    warning("Top left block is not square. Using full matrix.")
    C = full(A) \ B;
  elseif isscalar(B) && all(size(B) == 1)
    C = inv(A) * B;
  else
    [C1, C2] = blkmatrix_solve(A, B(1:size(A.A11,1),:), B(size(A.A11,1)+1:end,:));
    C = [C1; C2];
  end
end
end



