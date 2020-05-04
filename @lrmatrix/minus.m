function C = minus(A,B)

if any(size(A) ~= size(B))
  error("Incorrect dimensions for matrix substraction. Check that the dimensions of both matrices are matching.");
end

if isa(A, 'lrmatrix') && isa(B, 'lrmatrix')
  C = lrmatrix([A.U, -B.U], [A.V, B.V]);
elseif isa(A, 'lrmatrix')
  C = full(A) - B;
elseif isa(B, 'lrmatrix')
  C = A - full(B);
end

end