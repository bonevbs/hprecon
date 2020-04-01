function P = hprecon_fact(P, A)
  P = hprecon_symb_fact(P);
  [P, i] = hprecon_fact_rec(P, A, 0);
  assert(i == P.desc + 1);
end

function [P, i] = hprecon_fact_rec(P, A, i)
  if (isempty(P.Son1) && isempty(P.Son2))
    fprintf('factoring node %d\n', i)
    P = fact_leaf(P, A);
  elseif (~isempty(P.Son1) && ~isempty(P.Son2))
    [P.Son1, i] = hprecon_fact_rec(P.Son1, A, i);
    [P.Son2, i] = hprecon_fact_rec(P.Son2, A, i);
    fprintf('factoring node %d\n', i)
    P = fact_branch(P, A);
  else
    error('Elimination tree is not balanced. This is not supported.')
  end
   i = i + 1;
end

% assemble factorization from children
function P = fact_branch(P, A)
  % TODO: It might be completely unnecessary to permute the entries in the
  % outputted Schur complement. IT might be worth having a look at that


  % recover global indices in order to figure out the right entries for
  % off-diagonal blocks
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
    @(x) (x'*P.Son1.S)', @(x) A(idx2, idx1)'*x, ...
    @(x) A(idx1, idx2)'*x, @(x) (x'*P.Son2.S)', ...
    pi1, pb1, pi2, pb2, x);
  Aeval = @(i,j) blockeval( ...
    @(k,l) get(P.Son1.S, k, l), @(k,l) full(A(idx1(k), idx2(l))), ...
    @(k,l) full(A(idx2(k), idx1(l))), @(k,l) get(P.Son2.S, k, l), ...
    length(idx1), length(idx1), iperm(i) ,iperm(j));
  
  % compute the Schur complement through recompression
  cl = gen_cluster_rec([length(P.inter), length(P.inter) + length(P.bound)], hssoption('block-size'));
  hA = hss('handle', Afun, Afunt, Aeval, length(perm), length(perm), 'cluster', cl);
  
  % this is just some code to compare how well we are compressing
  Aref = [full(P.Son1.S), A(idx1, idx2); A(idx2, idx1), full(P.Son2.S)];
  Aref = Aref(iperm,iperm); norm(Aref - hA, 'fro')/norm(Aref, 'fro')
  % use the direct way to compute the reference matrix
  ii = setdiff(P.get_interior(), P.inter); bb = [idx1, idx2];
  Aref = A(bb,bb) - A(bb,ii) * (A(ii,ii) \ A(ii,bb));
  Aref = Aref(iperm,iperm); norm(Aref - hA, 'fro')/norm(Aref, 'fro')
  
  % only the top block and the left and right transforms need to be stored
  P.Aii = clean_hss(hA.A11);
  Abb = clean_hss(hA.A22);

  [P.RU,P.RV] = offdiag(hA, 'upper');
  [P.LU,P.LV] = offdiag(hA, 'lower');
  % compute intermediate matrix for efficient low-rank update
  K = (P.LV' * (P.Aii \ P.RU)) * P.RV';
  
  % compute the permutation to expose the finter block
  % TODO : Eventually we can get rid of this permutation step and also of
  % finter/fbound unless we develop a better way of assembling the
  % intermediate Schur complement which involves exposing the topmost block
  perm = []; iperm = [];
  perm = [P.finter, P.fbound]; iperm(perm) = 1:length(perm);
  
  % again using standard clusters, form the Schur complement in HSS form
  Sfun = @(x) Abb*x - P.LU*(K*x);
  Sfunt = @(x) Abb'*x - K'*(P.LU'*x);
  Seval = @(i,j) get(Abb, i, j) + P.LU(i,:)*K(:,j);
  P.S = hss('handle', @(x) Afun_perm(Sfun, iperm, x), ...
    @(x) Afun_perm(Sfun, iperm, x), ...
    @(i,j) Seval(perm(i), perm(j)), length(perm), length(perm));
  
  P.RU = P.Aii \ P.RU;
  P.LV = (P.LV' / P.Aii)';
end

% factor leafnodes directly
function P = fact_leaf(P, A)
  P.Aii = A(P.inter, P.inter);
  P.LU = A(P.bound, P.inter) / P.Aii;
  P.LV = speye(length(P.inter),length(P.inter));
  P.RV = (P.Aii \ A(P.inter, P.bound))';
  P.RU = speye(length(P.inter),length(P.inter));
  S      = A(P.bound, P.bound) - A(P.bound, P.inter) * (P.Aii \ A(P.inter, P.bound));
  cl     = gen_cluster_rec([length(P.finter), length(P.finter) + length(P.fbound)], hssoption('block-size'));
  % nonstandard blocking doesn't seem to be supported yet
  P.S    = hss(...
           [S(P.finter, P.finter), S(P.finter, P.fbound); S(P.fbound, P.finter), S(P.fbound, P.fbound)],...
           'cluster', cl);
  % for now only standard clusters
  %P.S    = hss([S(P.finter, P.finter), S(P.finter, P.fbound); S(P.fbound, P.finter), S(P.fbound, P.fbound)]);
end


% auxiliary functions to compute block matrix-vector multiplication
function z = blockfun_perm(A11fun, A12fun, A21fun, A22fun, pi1, pb1, pi2, pb2, x)
  z = zeros(size(x));
  z([pi1, pb1], :) = A11fun(x([pi1, pb1], :)) + A12fun(x([pi2, pb2], :));
  z([pi2, pb2], :) = A21fun(x([pi1, pb1], :)) + A22fun(x([pi2, pb2], :));
end

% auxiliary functions to compute entries in the blocked matrix
function e = blockeval(A11, A12, A21, A22, ni, nj, i ,j)
  e = zeros(length(i), length(j));
  is = (i <= ni); js = (j <= nj);
  if any(is) && any(js)
    e(is, js) = A11(i(is), j(js));
  end
  is = (i <= ni); js = (j > nj);
  if any(is) && any(js)
    e(is, js) = A12(i(is), j(js)-nj);
  end
  is = (i > ni); js = (j <= nj);
  if any(is) && any(js)
    e(is, js) = A21(i(is)-ni, j(js));
  end
  is = (i > ni); js = (j > nj);
  if any(is) && any(js)
    e(is, js) = A22(i(is)-ni, j(js)-nj);
  end
end

function x = Afun_perm(Afun, perm, x)
  x(perm, :) = Afun(x(perm,:));
end

% not necessary anymore
% wrapper function to sort entries before accessing them in the HSS matrix
%function out = sorted_access(A, i, j)
%  [~,is] = sort(i); [~,js] = sort(j);
%  out(is,js) = A(i(is), j(js));
%end