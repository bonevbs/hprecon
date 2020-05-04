function B = inv(A)
%INV Compute the inverse of an HSS matrix.

if size(A, 1) ~= size(A, 2)
    error('This matrix is not square');
end

B = A \ blkmatrix(eye(size(A.A11)),zeros(size(A.A12)),zeros(size(A.A21)),eye(size(A.A22)));

end
