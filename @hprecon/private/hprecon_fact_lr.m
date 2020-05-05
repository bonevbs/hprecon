function hprecon_fact_lr(P, A)
  P = hprecon_symb_fact(P);
  [P, i] = hprecon_fact_rec(P, A, 1, 1, hpreconoption('levels'));
  %assert(i == P.desc + 1);
end

function [P, i] = hprecon_fact_rec(P, A, i, l, lmax)
  if l == lmax || (isempty(P.Son1) && isempty(P.Son2))
    %fprintf('HPRECON_FACT: factoring node %d\n', i)
    P = fact_leaf(P, A);
  elseif (l < lmax || lmax == -1) && (~isempty(P.Son1) && ~isempty(P.Son2))
    [P.Son1, i] = hprecon_fact_rec(P.Son1, A, i, l+1, lmax);
    [P.Son2, i] = hprecon_fact_rec(P.Son2, A, i, l+1, lmax);
    %fprintf('HPRECON_FACT: factoring node %d\n', i)
    P = fact_branch(P, A);
  else
    error('HPRECON_FACT: Something went wrong. This is not supported.')
  end
   i = i + 1;
end

% assemble factorization from children
function P = fact_branch(P, A)

  % recover global indices in order to figure out the right entries for
  % off-diagonal blocks
  inter1 = P.Son1.bound(P.Son1.finter); 
  inter2 = P.Son2.bound(P.Son2.finter);
  bound1 = P.Son1.bound(P.Son1.fbound);
  bound2 = P.Son2.bound(P.Son2.fbound);
  
  ni1 = length(inter1); ni2 = length(inter2);
  nb1 = length(bound1); nb2 = length(bound2);
  
  % datastructure to form the blockinverse
  P.Aii = blkmatrix(clean_hss(P.Son1.S.A11), A(inter1,inter2), A(inter2,inter1), clean_hss(P.Son2.S.A11));
  Abb = blkmatrix(clean_hss(P.Son1.S.A22), A(bound1,bound2), A(bound2,bound1), clean_hss(P.Son2.S.A22));
  
  %Aii = [full(P.Aii.A11), full(P.Aii.A12); full(P.Aii.A21), full(P.Aii.A22)];
  %Abb = [full(clean_hss(P.Son1.S.A22)), full(A(bound1,bound2)); full(A(bound2,bound1)), full(clean_hss(P.Son2.S.A22))];

  % create function handle to apply the left Gauss transform
  Abi11 = lrmatrix(); [Abi11.U, Abi11.V] = offdiag(P.Son1.S,'lower');
  Abi12 = A(bound1,inter2);
  Abi21 = A(bound2,inter1);
  Abi22 = lrmatrix(); [Abi22.U, Abi22.V] = offdiag(P.Son2.S,'lower');
  Abi = blkmatrix(Abi11, Abi12, Abi21, Abi22);
  [U,S,V] = svd_from_random_sampling(@(x) Abi*(P.Aii\x), @(x) ((x'*Abi)/P.Aii)', ni1+ni2, 1e-9);
  P.L = lrmatrix(U*S, V);
  
  % apply the inverse to Aib to form the right Gauss transform
  Aib11 = lrmatrix(); [Aib11.U, Aib11.V] = offdiag(P.Son1.S,'upper');
  Aib12 = A(inter1,bound2);
  Aib21 = A(inter2,bound1);
  Aib22 = lrmatrix(); [Aib22.U, Aib22.V] = offdiag(P.Son2.S,'upper');
  Aib = blkmatrix(Aib11, Aib12, Aib21, Aib22);
  [U,S,V] = svd_from_random_sampling(@(x) P.Aii\(Aib*x), @(x) ((x'/P.Aii)*Aib)', nb1+nb2, 1e-9);
  P.R = lrmatrix(U*S, V);
  
  K = Abi * P.R;

  % the dense version
  S = full(Abb) - K;
  S = S([P.finter,P.fbound],[P.finter,P.fbound]);
  %P.S = hss(P.S);
  %P.S = hss('handle', @(x) P.S*x, @(x) P.S'*x, @(i,j) P.S(i,j), size(P.S,1), size(P.S,2));
  %P.S = hss(S);
  P.S = hss(S,'cluster', gen_cluster_rec([length(P.finter),length(P.finter)+length(P.fbound)],hssoption('block-size')));

  %fprintf("HPRECON_FACT_HSS: hss rank of the Schur complement: %d\n", hssrank(P.S));
end

% factor leafnodes directly
function P = fact_leaf(P, A)
  inter = P.get_interior();
  bound = P.bound;
  P.Aii = A(inter, inter);
  P.L = A(bound, inter) / P.Aii;
  P.R = P.Aii \ A(inter, bound);
  S = A(bound, bound) - A(bound, inter) * (P.Aii \ A(inter, bound));
  S = S([P.finter,P.fbound],[P.finter,P.fbound]);
  %P.S = hss(S);
  P.S = hss(S,'cluster', gen_cluster_rec([length(P.finter),length(P.finter)+length(P.fbound)],hssoption('block-size')));
end
