% this is the brute force version of the factorization, utilizing full
% matrices to assemble the entire thing

function maxk = hprecon_fact(P, A, maxk)
  P = hprecon_symb_fact(P);
  maxk = 0;
  [P, i, maxk] = hprecon_fact_rec(P, A, 0, maxk);
  assert(i == P.desc + 1);
end

function [P, i, maxk] = hprecon_fact_rec(P, A, i, maxk)
  if (isempty(P.Son1) && isempty(P.Son2))
    %fprintf('HPRECON_FACT: factoring node %d\n', i)
    [P, maxk] = fact_leaf(P, A, maxk);
  elseif (~isempty(P.Son1) && ~isempty(P.Son2))
    [P.Son1, i, maxk] = hprecon_fact_rec(P.Son1, A, i, maxk);
    [P.Son2, i, maxk] = hprecon_fact_rec(P.Son2, A, i, maxk);
    %fprintf('HPRECON_FACT: factoring node %d\n', i)
    [P, maxk] = fact_branch(P, A, maxk);
  else
    error('HPRECON_FACT: Elimination tree is not balanced. This is not supported.')
  end
   i = i + 1;
end

% assemble factorization from children
function [P, maxk] = fact_branch(P, A, maxk)
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
  
  P.Aii = struct;
  P.Aii.S11 = hss(get(P.Son1.S, P.Son1.finter, P.Son1.finter));
  P.Aii.A12 = A(inter1,inter2);
  P.Aii.A21 = A(inter2,inter1);
  P.Aii.S22 = hss(get(P.Son2.S, P.Son2.finter, P.Son2.finter));

  function x = Lfun(x)
    [x1, x2] = blocksolve_right(P.blAii.S11, P.blAii.A12, P.blAii.A21, P.blAii.S22, x(1:length(inter1)), x(end-length(inter1)+1:end));
  end
  
  
  % only the top block and the left and right transforms need to be stored
  %P.invAii = hodlr(inv(Aii));
  %hodlrrank(P.invAii)
  %P.invAii = inv(Aii);
  %[U,S,V] = tsvd(full(Abi/Aii), 1e-6);
  [U,S,V] = svd_from_random_sampling(@(x) full(Abi/Aii)*x, @(x) full(Abi/Aii)'*x, length(P.inter), 1e-4);
  P.LU = U*S; P.LV = V;
  norm(full(Abi/Aii) - P.LU * P.LV', 'fro')/norm(full(Abi/Aii), 'fro');
  %[U,S,V] = tsvd(full(Aii\Aib), 1e-6);
  [U,S,V] = svd_from_random_sampling(@(x) full(Aii\Aib)*x, @(x) full(Aii\Aib)'*x, length(P.bound), 1e-4);
  P.RU = U*S; P.RV = V;

%   P.Aii = Aii;
%   P.LU = Abi / P.Aii;
%   P.LV = speye(length(P.inter),length(P.inter));
%   P.RV = (P.Aii \ Aib)';
%   P.RU = speye(length(P.inter),length(P.inter));
  
  

  %P.RU = hA.U12; P.RV = hA.V12;
  %P.LU = hA.U21; P.LV = hA.V21;
  % compute intermediate matrix for efficient low-rank update
  K = P.LU * (P.LV' * (Aii * P.RU)) * P.RV';
  
  
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
  %fprintf("HPRECON_FACT: hss rank of the Schur complement: %d\n", hssrank(P.S));
  if (size(P.LU,2) > maxk) && (size(P.RU,2) > maxk) && (hssrank(P.S) > maxk)
    maxk = max([size(P.LU,2), size(P.RU,2), hssrank(P.S)]);
  end
  
  
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
function [P, maxk] = fact_leaf(P, A, maxk)
  %P.invAii = inv(A(P.inter, P.inter));
  P.Aii = A(P.inter, P.inter);
  P.LU = A(P.bound, P.inter) / P.Aii;
  P.LV = speye(length(P.inter),length(P.inter));
  P.RV = (P.Aii \ A(P.inter, P.bound))';
  P.RU = speye(length(P.inter),length(P.inter));
  P.S      = hss(A(P.bound, P.bound) - A(P.bound, P.inter) * (P.Aii \ A(P.inter, P.bound)));
  %P.S      = A(P.bound, P.bound) - A(P.bound, P.inter) * (P.Aii \ A(P.inter, P.bound));
  if (hssrank(P.S) > maxk)
    maxk = hssrank(P.S);
  end
end