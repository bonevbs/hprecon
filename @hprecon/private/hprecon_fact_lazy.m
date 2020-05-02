function hprecon_fact_full(P, A)
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
  
  Aii = [get(P.Son1.S, 1:length(inter1), 1:length(inter1)), A(inter1,inter2); A(inter2,inter1), get(P.Son2.S, 1:length(inter2), 1:length(inter2))];
  Aib = [get(P.Son1.S, 1:length(inter1), length(inter1)+1:length(inter1)+length(bound1)), A(inter1,bound2); A(inter2,bound1), get(P.Son2.S, 1:length(inter2), length(inter2)+1:length(inter2)+length(bound2))];
  Abi = [get(P.Son1.S, length(inter1)+1:length(inter1)+length(bound1), 1:length(inter1)), A(bound1,inter2); A(bound2,inter1), get(P.Son2.S, length(inter2)+1:length(inter2)+length(bound2), 1:length(inter2))];
  Abb = [get(P.Son1.S, length(inter1)+1:length(inter1)+length(bound1), length(inter1)+1:length(inter1)+length(bound1)), A(bound1,bound2); A(bound2,bound1), get(P.Son2.S, length(inter2)+1:length(inter2)+length(bound2), length(inter2)+1:length(inter2)+length(bound2))];

  % datastructure to form the blockinverse
  P.Aii = struct;
  P.Aii.A11 = clean_hss(P.Son1.S.A11);
  P.Aii.A12 = A(inter1,inter2);
  P.Aii.A21 = A(inter2,inter1);
  P.Aii.A22 = clean_hss(P.Son2.S.A11);

  [invAii11, invAii21] = blksolve_right(P.Aii, eye(size(P.Aii.A11)), zeros(size(P.Aii.A21)));
  [invAii12, invAii22] = blksolve_right(P.Aii, zeros(size(P.Aii.A12)), eye(size(P.Aii.A22)));
  
  P.invAii = [invAii11, invAii12; invAii21, invAii22];   
  P.LU = Abi * P.invAii;
  P.LV = speye(length(P.inter),length(P.inter));
  P.RV = (P.invAii * Aib)';
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

% not necessary anymore
% wrapper function to sort entries before accessing them in the HSS matrix
%function out = sorted_access(A, i, j)
%  [~,is] = sort(i); [~,js] = sort(j);
%  out(is,js) = A(i(is), j(js));
%end