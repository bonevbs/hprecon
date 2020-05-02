function C = mtimes(A,B)

if size(A,2) ~= size(B,1)
  error("Incorrect dimensions for matrix multiplication. Check that the number of columns in the first matrix matches the number of rows in the second matrix.");
end

if isa(A, 'lrmatrix') && isa(B, 'lrmatrix')
  C = lrmatrix(A.U, B.V*(B.U'*A.V));
elseif isa(A, 'lrmatrix')
  C = lrmatrix(A.U, (A.V'*B)');
elseif isa(B, 'lrmatrix')
  C = lrmatrix(A*B.U, B.V);
end

end