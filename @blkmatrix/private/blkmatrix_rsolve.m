function [x1, x2] = blkmatrix_rsolve(A, b1, b2)
x1 = b1/A.A11;
x2 = b2 - x1*A.A12;
if isempty(A.S22)
  if isa(A.A11, 'hss') && isa(A.A12, 'hss') && isa(A.A21, 'hss') && isa(A.A22, 'hss')
    A.S22 = A.A22 - hss(A.A21*(full(A.A11)\A.A12), 'cluster', cluster(A.A22));
    assert(isa(A.S22, 'hss'), "Something went wrong. Schur complement is not HSS.")
  elseif isa(A.A11, 'hss') && isa(A.A22, 'hss')
    % this should be replaced by a method that uses randomized sampling
    A.S22 = hss(A.A22 - A.A21*(A.A11\A.A12));
    % this is the lazy implementation of the direct sampling method
    %S = A.A22 - A.A21*(A.A11\A.A12);
    %P.S = build_hss_tree(size(S,1), size(S,2), hssoption('block-size'), [], []);
    %A.S22 = hss_martinsson_adaptive(P.S, @(x) S*x, @(x) S'*x, @(i,j) S(i,j), size(S,1), size(S,2), 10);
  else
    A.S22 = A.A22 - A.A21*(A.A11\A.A12);
  end
end
x2 = x2/A.S22;
x1 = x1 - (x2*A.A21)/A.A11;
end