function [P, i] = compute_factorization_rec(P, A, i)
  if (isempty(P.Son1) && isempty(P.Son2))
    fprintf('factoring node %d\n', i)
    P = factor_leafnode(P, A);
  elseif (~isempty(P.Son1) && ~isempty(P.Son2))
    [P.Son1, i] = compute_factorization_rec(P.Son1, A, i);
    [P.Son2, i] = compute_factorization_rec(P.Son2, A, i);
    fprintf('factoring node %d\n', i)
    P = factor_branchnode(P, A);
  else
    error('Elimination tree is not balanced. This is not supported.')
  end
   i = i + 1;
end

% assemble factorization from children
function P = factor_branchnode(P, A)

  inter1 = P.Son1.bound(P.Son1.finter);
  inter2 = P.Son2.bound(P.Son2.finter);
  bound1 = P.Son1.bound(P.Son1.fbound);
  bound2 = P.Son2.bound(P.Son2.fbound);
  idx1 = [inter1, bound1];
  idx2 = [inter2, bound2];

  % first step - assemble & compress the combined Schur complement
  ni1 = length(P.Son1.finter);
  ni2 = ni1 + length(P.Son2.finter);
  nb1 = ni2 + length(P.Son1.fbound);
  nb2 = nb1 + length(P.Son2.fbound);
  pi1 = 1:ni1; pi2 = ni1+1:ni2; pb1 = ni2+1:nb1; pb2 = nb1+1:nb2;
  perm = [pi1, pb1, pi2, pb2]; iperm(perm) = 1:length(perm);
  
  % the interface is [pi1, pi2, pb1, pb2]
  % the hidden matrix is structured perm = [pi1, pb1, pi2, pb2]
  % x(perm) allows us to convert a vector from the interface format to the
  % internal one.
  Afun = @(x) blockfun_perm( ...
    @(x) P.Son1.S*x, @(x) A(idx1, idx2)*x, ...
    @(x) A(idx2, idx1)*x, @(x) P.Son2.S*x, ...
    pi1, pb1, pi2, pb2, x);
  Afunt = @(x) blockfun_perm( ...
    @(x) (x'*P.Son1.S)', @(x) A(idx2,idx1)'*x, ...
    @(x) A(idx1, idx2)'*x, @(x) (x'*P.Son2.S)', ...
    pi1, pb1, pi2, pb2, x);
  
  
  
  
  Aeval = @(i,j) blockeval( ...
    @(k,l) full(P.Son1.S(k,l)), @(k,l) full(A(idx1(k),idx2(l))), ...
    @(k,l) full(A(idx2(k), idx1(l))), @(k,l) full(P.Son2.S(k,l)), ...
    length(idx1), length(idx1), perm(i) ,perm(j));
  % the output vector should be orderered [inter]
  
  % compute the Schur complement through recompression
  cl = gen_cluster_rec([length(P.inter), length(P.inter) + length(P.bound)], hssoption('block-size'));
  A = hss('handle', Afun, Afunt, Aeval, length(perm), length(perm));
end

% factor leafnodes directly
function P = factor_leafnode(P, A)
  P.Aii  = hss(inv(A(P.inter, P.inter)));
  %P.invAB = P.invA * A(P.inter, P.bound);
  %P.CinvA = A(P.bound, P.inter) * P.invA;p
  %P.S     = A(P.bound, P.bound) - A(P.bound, P.inter) * P.invAB;
  S      = A(P.bound, P.bound) - A(P.bound, P.inter) * (P.Aii \ A(P.inter, P.bound));
  cl     = gen_cluster_rec([length(P.finter), length(P.finter) + length(P.fbound)], hssoption('block-size'));
  P.S    = hss(...
           [S(P.finter, P.finter), S(P.finter, P.fbound); S(P.fbound, P.finter), S(P.fbound, P.fbound);], ...
           'cluster', cl);
end


% auxiliary functions to compute block matrix-vector multiplication
function z = blockfun_perm(A11fun, A12fun, A21fun, A22fun, pi1, pb1, pi2, pb2, x)
  z = zeros(size(x));
  z([pi1, pb1], :) = A11fun(x([pi1, pb1], :)) + A12fun(x([pi2, pb2], :));
  z([pi2, pb2], :) = A21fun(x([pi1, pb1], :)) + A22fun(x([pi2, pb2], :));
end

% auxiliary functions to compute entries in the blocked matrix
function e = blockeval(A11eval, A12eval, A21eval, A22eval, ni, nj, i ,j)
  e = zeros(length(i), length(j));
  e((i <= ni), (j <= nj)) = sorted_access(A11eval, i(i <= ni), j(j <= nj));
  e((i <= ni), (j > nj) ) = A12eval(i(i <= ni), j(j > nj)-nj);
  e( (i > ni), (j <= nj)) = A21eval(i(i > ni)-ni, j(j <= nj));
  e( (i > ni), (j > nj) ) = sorted_access(A22eval, i(i > ni)-ni, j(j > nj)-nj );
end

% wrapper function to sort entries before accessing them in the HSS matrix
function out = sorted_access(Aeval, i, j)
  [~,is] = sort(i); [~,js] = sort(j);
  out(is,js) = Aeval(i(is), j(js));
end

% this is not efficient in theory as it scales as O(N). As we need to
% access O(N) entries, this is what makes the algorithm quadratic as of now
%function out = permuted_access(Afun, i, j)
%end