% this is the brute force version of the factorization, utilizing full
% matrices to assemble the entire thing

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
  
  Aii = [get(P.Son1.S, P.Son1.finter, P.Son1.finter), A(inter1,inter2); A(inter2,inter1), get(P.Son2.S, P.Son2.finter, P.Son2.finter)];
  Aib = [get(P.Son1.S, P.Son1.finter, P.Son1.fbound), A(inter1,bound2); A(inter2,bound1), get(P.Son2.S, P.Son2.finter, P.Son2.fbound)];
  Abi = [get(P.Son1.S, P.Son1.fbound, P.Son1.finter), A(bound1,inter2); A(bound2,inter1), get(P.Son2.S, P.Son2.fbound, P.Son2.finter)];
  Abb = [get(P.Son1.S, P.Son1.fbound, P.Son1.fbound), A(bound1,bound2); A(bound2,bound1), get(P.Son2.S, P.Son2.fbound, P.Son2.fbound)];
  
%   Aii = [P.Son1.S(P.Son1.finter, P.Son1.finter), A(inter1,inter2); A(inter2,inter1), P.Son2.S(P.Son2.finter, P.Son2.finter)];
%   Aib = [P.Son1.S(P.Son1.finter, P.Son1.fbound), A(inter1,bound2); A(inter2,bound1), P.Son2.S(P.Son2.finter, P.Son2.fbound)];
%   Abi = [P.Son1.S(P.Son1.fbound, P.Son1.finter), A(bound1,inter2); A(bound2,inter1), P.Son2.S(P.Son2.fbound, P.Son2.finter)];
%   Abb = [P.Son1.S(P.Son1.fbound, P.Son1.fbound), A(bound1,bound2); A(bound2,bound1), P.Son2.S(P.Son2.fbound, P.Son2.fbound)];

  %Aref = [Aii, Aib; Abi, Abb];
  %cl = gen_cluster_rec([length(inter1), length(inter1) + length(inter2)], hssoption('block-size'));
  %hAii = hodlr(Aii, 'cluster', cl);

  % use the direct way to compute the reference matrix
  %ii = setdiff(P.get_interior(), P.inter); bb = [inter1, inter2, bound1, bound2];
  %Aref = A(bb,bb) - A(bb,ii) * (A(ii,ii) \ A(ii,bb));
  %norm(Aref - hA, 'fro')/norm(Aref, 'fro')
  
  
  
  % only the top block and the left and right transforms need to be stored
  %P.invAii = hodlr(inv(Aii));
  %hodlrrank(P.invAii)
  P.Aii = Aii;
  P.invAii = inv(Aii);
  %[U,S,V] = tsvd(full(Abi) / P.Aii, 1e-3);
  [U,S,V] = svd_from_random_sampling(@(x) full(Abi/P.Aii)*x, @(x) full(Abi/P.Aii)'*x, length(P.inter), 1e-4);
  P.LU = U*S; P.LV = V;
  norm(full(Abi/P.Aii) - P.LU * P.LV', 'fro')/norm(full(Abi/P.Aii), 'fro');
  %[U,S,V] = tsvd(full(P.Aii\Aib), 1e-3);
  [U,S,V] = svd_from_random_sampling(@(x) full(P.Aii\Aib)*x, @(x) full(P.Aii\Aib)'*x, length(P.bound), 1e-4);
  
  P.RU = U*S; P.RV = V;

%   P.Aii = Aii;
%   P.LU = Abi / P.Aii;
%   P.LV = speye(length(P.inter),length(P.inter));
%   P.RV = (P.Aii \ Aib)';
%   P.RU = speye(length(P.inter),length(P.inter));
  
  

  %P.RU = hA.U12; P.RV = hA.V12;
  %P.LU = hA.U21; P.LV = hA.V21;
  % compute intermediate matrix for efficient low-rank update
  K = P.LU * (P.LV' * (P.Aii * P.RU)) * P.RV';
  
  
%   % again using standard clusters, form the Schur complement in HSS form
%   Sfun = @(x) Abb*x - P.LU*(K*x);
%   Sfunt = @(x) Abb'*x - K'*(P.LU'*x);
%   Seval = @(i,j) get(Abb, i, j) + P.LU(i,:)*K(:,j);
%   P.S = hss('handle', @(x) Afun_perm(Sfun, iperm, x), ...
%     @(x) Afun_perm(Sfun, iperm, x), ...
%     @(i,j) Seval(perm(i), perm(j)), length(perm), length(perm));
  
  % the dense version
  P.S = Abb - K;
  P.S = hss(P.S);
  
  %P.RU = P.Aii \ P.RU;
  %P.LV = (P.LV' / P.Aii)';
  hssrank(P.S)
  
  % overwrite everything with the dense version to make sure everything
  % else works correctly
%   P.Aii = Aref([pi1,pi2], [pi1,pi2]);
%   P.LU = A([pb1,pb2], [pi1,pi2]) / P.Aii;
%   P.LV = speye(length(P.inter),length(P.inter));
%   P.RV = (P.Aii \ A([pi1,pi2], [pb1,pb2]))';
%   P.RU = speye(length(P.inter),length(P.inter));
%   P.S = Aref([pb1, pb2], [pb1, pb2]) - Aref([pb1, pb2], [pi1, pi2]) *(Aref([pi1, pi2], [pi1, pi2]) \ Aref([pi1, pi2], [pb1, pb2]));
%   P.S = P.S([P.finter, P.fbound], [P.finter, P.fbound]);
end

% factor leafnodes directly
function P = fact_leaf(P, A)
  P.invAii = inv(A(P.inter, P.inter));
  P.Aii = A(P.inter, P.inter);
  P.LU = A(P.bound, P.inter) / P.Aii;
  P.LV = speye(length(P.inter),length(P.inter));
  P.RV = (P.Aii \ A(P.inter, P.bound))';
  P.RU = speye(length(P.inter),length(P.inter));
  P.S      = hss(A(P.bound, P.bound) - A(P.bound, P.inter) * (P.Aii \ A(P.inter, P.bound)));
  %P.S      = A(P.bound, P.bound) - A(P.bound, P.inter) * (P.Aii \ A(P.inter, P.bound));
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