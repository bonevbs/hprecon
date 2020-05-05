function C = mtimes(A,B)

if size(A,2) ~= size(B,1)
  error("Incorrect dimensions for matrix multiplication. Check that the number of columns in the first matrix matches the number of rows in the second matrix.");
end

if isa(A, 'blkmatrix') && isa(B, 'blkmatrix')
  if size(A.A11,2) == size(B.A11,1)
    C = blkmatrix(A.A11*B.A11 + A.A12*B.A21, A.A11*B.A12 + A.A12*B.A22, ...
      A.A21*B.A11 + A.A22*B.A21, A.A21*B.A12 + A.A22*B.A22);
  else
    warning("Matrix blocking doesn't align. This may be inefficient")
    C = A * full(B);
  end
elseif isa(A, 'blkmatrix')
  C = [A.A11*B(1:size(A.A11,2), :) + A.A12*B(size(A.A11,2)+1:end, :);...
    A.A21*B(1:size(A.A11,2), :) + A.A22*B(size(A.A11,2)+1:end, :)];
elseif isa(B, 'blkmatrix')
  C = [A(:,1:size(B.A11,1))*B.A11 + A(:,size(B.A11,1)+1:end)*B.A21,...
    A(:,1:size(B.A11,1))*B.A12 + A(:,size(B.A11,1)+1:end)*B.A22];
end

end