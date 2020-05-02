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
  P.Aii = struct;
  P.Aii.A11 = clean_hss(P.Son1.S.A11);
  P.Aii.A12 = A(inter1,inter2);
  P.Aii.A21 = A(inter2,inter1);
  P.Aii.A22 = clean_hss(P.Son2.S.A11);
  
  Aii = [full(P.Aii.A11), full(P.Aii.A12); full(P.Aii.A21), full(P.Aii.A22)];
  Abb = [full(clean_hss(P.Son1.S.A22)), full(A(bound1,bound2)); full(A(bound2,bound1)), full(clean_hss(P.Son2.S.A22))];

  
  % create function handle to apply the left Gauss transform
  [U1, V1] = offdiag(P.Son1.S,'lower'); [U2, V2] = offdiag(P.Son2.S,'lower');
  L12 = A(bound1,inter2); L21 = A(bound2,inter1);
  function x = Lfun(x)
    [x1, x2] = blksolve_right(P.Aii, x(1:ni1,:), x(ni1+1:ni1+ni2,:));
    x = [U1*(V1'*x1) + L12*x2; L21*x1 + U2*(V2'*x2)];
  end
  function x = Lfunt(x)
    x = x';
    x1 = (x(:,1:nb1)*U1)*V1' + x(:,nb1+1:nb1+nb2)*L21;
    x2 = x(:,1:nb1)*L12 + (x(:,nb1+1:nb1+nb2)*U2)*V2';
    [x1, x2] = blksolve_left(P.Aii, x1, x2);
    x = [x1, x2]; x = x';
  end
  [U,S,V] = svd_from_random_sampling(@(x) Lfun(x), @(x) Lfunt(x), ni1+ni2, 1e-9);
  P.LU = U*S; P.LV = V;
  
  % apply the inverse to Abi to form left Gauss transform
%   [U, V] = offdiag(P.Son1.S,'lower'); L11 = U*V'; L12 = A(bound1,inter2);
%   [U, V] = offdiag(P.Son2.S,'lower'); L22 = U*V'; L21 = A(bound2,inter1);
%   [L11, L12] = blksolve_left(P.Aii, L11, L12);
%   [L21, L22] = blksolve_left(P.Aii, L21, L22);
%   P.LU = [L11,L12;L21,L22];
%   P.LV = speye(length(P.inter),length(P.inter));
  
  % apply the inverse to Aib to form the right Gauss transform
  [U, V] = offdiag(P.Son1.S,'upper'); R11 = U*V'; R21 = A(inter2,bound1);
  [U, V] = offdiag(P.Son2.S,'upper'); R22 = U*V'; R12 = A(inter1,bound2);
  [R11, R21] = blksolve_right(P.Aii, R11, R21);
  [R12, R22] = blksolve_right(P.Aii, R12, R22);
  P.RV = [R11,R12;R21,R22]';
  P.RU = speye(length(P.inter),length(P.inter));
  
  K = P.LU * (P.LV' * (Aii * P.RU)) * P.RV';

  % the dense version
  S = Abb - K;
  S = S([P.finter,P.fbound],[P.finter,P.fbound]);
  %P.S = hss(P.S);
  %P.S = hss('handle', @(x) P.S*x, @(x) P.S'*x, @(i,j) P.S(i,j), size(P.S,1), size(P.S,2));
  %P.S = hss(S);
  P.S = hss(S,'cluster', gen_cluster_rec([length(P.finter),length(P.finter)+length(P.fbound)],hssoption('block-size')));
  P.invAii = [];

  %fprintf("HPRECON_FACT_HSS: hss rank of the Schur complement: %d\n", hssrank(P.S));
end

% factor leafnodes directly
function P = fact_leaf(P, A)
  inter = P.get_interior();
  bound = P.bound;
  P.Aii = A(inter, inter);
  P.LU = A(bound, inter) / P.Aii;
  P.LV = speye(length(inter),length(inter));
  P.RV = (P.Aii \ A(inter, bound))';
  P.RU = speye(length(inter),length(inter));
  S  = A(bound, bound) - A(bound, inter) * (P.Aii \ A(inter, bound));
  S = S([P.finter,P.fbound],[P.finter,P.fbound]);
  %P.S = hss(S);
  P.S = hss(S,'cluster', gen_cluster_rec([length(P.finter),length(P.finter)+length(P.fbound)],hssoption('block-size')));
end
