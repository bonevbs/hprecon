function P = hprecon_fact_hss(P, A)
  P = hprecon_symb_fact(P);
  if hpreconoption('levels') < 0
    L = maxlevel(P) + hpreconoption('levels');
  else
    L = hpreconoption('levels');
  end
  [P, i] = hprecon_fact_rec(P, A, 1, 1, L);
  %assert(i == P.desc + 1);
end

function [P, i] = hprecon_fact_rec(P, A, i, l, lmax)
  if l == lmax || (isempty(P.Son1) && isempty(P.Son2))
    fprintf('HPRECON_FACT: factoring node %d\n', i)
    P = fact_leaf(P, A);
  elseif (l < lmax || lmax == 0) && (~isempty(P.Son1) && ~isempty(P.Son2))
    [P.Son1, i] = hprecon_fact_rec(P.Son1, A, i, l+1, lmax);
    [P.Son2, i] = hprecon_fact_rec(P.Son2, A, i, l+1, lmax);
    fprintf('HPRECON_FACT: factoring node %d\n', i)
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
  
  % blkmatrix form of Aii
  P.Aii = blkmatrix();
  P.Aii.A11 = clean_hss(P.Son1.S.A11);
  P.Aii.A22 = clean_hss(P.Son2.S.A11);
  [r1, c1] = cluster(P.Aii.A11);
  [r2, c2] = cluster(P.Aii.A22);
  P.Aii.A12 = hss(A(inter1,inter2), 'cluster', r1, c2);
  P.Aii.A21 = hss(A(inter2,inter1), 'cluster', r2, c1);
  %max(hssrank(P.Aii.A12), hssrank(P.Aii.A21))
  %norm(P.Aii.A12 - A(inter1,inter2), 'fro')/norm(A(inter1,inter2), 'fro')
  %norm(P.Aii.A21 - A(inter2,inter1), 'fro')/norm(A(inter2,inter1), 'fro')

  %P.Aii = blkmatrix(clean_hss(P.Son1.S.A11), A(inter1,inter2), A(inter2,inter1), clean_hss(P.Son2.S.A11));
 
  % Abb
  Abb = blkmatrix(clean_hss(P.Son1.S.A22), A(bound1,bound2), A(bound2,bound1), clean_hss(P.Son2.S.A22));
  
  % create function handle to apply the left Gauss transform
  Abi11 = lrmatrix(); [Abi11.U, Abi11.V] = offdiag(P.Son1.S,'lower');
  Abi12 = A(bound1,inter2);
  Abi21 = A(bound2,inter1);
  Abi22 = lrmatrix(); [Abi22.U, Abi22.V] = offdiag(P.Son2.S,'lower');
  Abi = blkmatrix(Abi11, Abi12, Abi21, Abi22);
  
  % apply the inverse to Aib to form the right Gauss transform
  Aib11 = lrmatrix(); [Aib11.U, Aib11.V] = offdiag(P.Son1.S,'upper');
  Aib12 = A(inter1,bound2);
  Aib21 = A(inter2,bound1);
  Aib22 = lrmatrix(); [Aib22.U, Aib22.V] = offdiag(P.Son2.S,'upper');
  Aib = blkmatrix(Aib11, Aib12, Aib21, Aib22);
  
  % form and left Gauss transforms
  P.L = lrmatrix(@(x) Abi*(P.Aii\x), @(x) ((x'*Abi)/P.Aii)', [nb1+nb2, ni1+ni2]);
  P.R = lrmatrix(@(x) P.Aii\(Aib*x), @(x) ((x'/P.Aii)*Aib)', [ni1+ni2, nb1+nb2]);
  
  K = Abi * P.R;
  perm = [P.finter,P.fbound];
  iperm(perm) = 1:length(perm);
  
  % these functions are necessary for 
  function x = Sfun(x)
    x(iperm,:) = Abb * x(iperm,:) - K * x(iperm,:);
  end
  function x = Sfunt(x)
    x(iperm,:) = (x(iperm,:)' * Abb - x(iperm,:)' * K)';
  end
  function x = Seval(i,j)
    x = Abb(perm(i),perm(j)) - K(perm(i),perm(j));
  end
  
  if strcmp(hpreconoption('merging-algorithm'), 'martinsson')
    cl = gen_cluster_rec([length(P.finter),length(P.finter)+length(P.fbound)],hssoption('block-size'));
    %P.S = hss('handle', @Sfun, @Sfunt, @Seval, size(Abb,1), size(Abb,2),...
    %  'cluster', gen_cluster_rec([length(P.finter),length(P.finter)+length(P.fbound)],hssoption('block-size')));
    S(iperm,iperm) = full(Abb) - full(K);
    %P.S = hss('handle', @(x) S*x, @(x) S'*x, @(i,j) S(i,j), size(Abb,1), size(Abb,2),'cluster', cl);
    P.S = build_hss_tree(size(S,1), size(S,2), hssoption('block-size'), cl, cl);
    P.S = hss_martinsson_adaptive(P.S, @Sfun, @Sfunt, @Seval, size(Abb,1), size(Abb,2), 10);
  elseif strcmp(hpreconoption('merging-algorithm'), 'low-rank')
    %B = lrmatrix([sparse(nb1,nb1), A(bound1,bound2), A(bound2,bound1), sparse(nb2,nb2)]) + K;
    %B = hss('low-rank', B.U, B.V);
    error('Unsupported as of now')
  elseif strcmp(hpreconoption('merging-algorithm'), 'direct')
    S(iperm,iperm) = full(Abb) - full(K);
    P.S = hss(S,'cluster', gen_cluster_rec([length(P.finter),length(P.finter)+length(P.fbound)],hssoption('block-size')));
    %P.S = proper(P.S);
    %P.S = compress(P.S);
  else
    error('Unknown merging algorithm.')
  end

  %S = S([P.finter,P.fbound],[P.finter,P.fbound]);
  %P.S = hss('handle', @(x) P.S*x, @(x) P.S'*x, @(i,j) P.S(i,j), size(P.S,1), size(P.S,2));

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
  P.S = hss(S,'cluster', gen_cluster_rec([length(P.finter),length(P.finter)+length(P.fbound)],hssoption('block-size')));
  %P.S = compress(P.S);
end
