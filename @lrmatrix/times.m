function C = times(A, B)

if (isfloat(A) && isscalar(A)) || (isfloat(B) && isscalar(B))
    C = A * B;
elseif isa(A, 'lrmatrix')
    C = full(A) .* B;
elseif isa(B, 'lrmatrix')
    C = A .* full(B);
else
    error('A .* B: Unsupported case');
end
end