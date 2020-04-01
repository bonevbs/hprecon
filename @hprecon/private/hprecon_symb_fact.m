function [P] = hprecon_symb_fact(P)
  [P, i] = hprecon_symb_fact_rec(P, 0);
  assert(i == P.desc + 1);
  % fix finter and fbound for the top node
  P.finter = 1:length(P.bound);
end

function [P, i] = hprecon_symb_fact_rec(P, i)
  if (isempty(P.Son1) && isempty(P.Son2))
    % nada
  elseif (~isempty(P.Son1) && ~isempty(P.Son2))
    [P.Son1, i] = hprecon_symb_fact_rec(P.Son1, i);
    [P.Son2, i] = hprecon_symb_fact_rec(P.Son2, i);
    %fprintf('factoring node %d\n', i)
    P = symbfact(P);
  else
    error('Elimination tree is not balanced. This is not supported.')
  end
  i = i + 1;
end

% assemble factorization from children
function P = symbfact(P)
  % re-organize indices
  inter1 = setdiff(P.Son1.bound, P.bound, 'stable');
  inter2 = setdiff(P.Son2.bound, P.bound, 'stable');
  bound1 = setdiff(P.Son1.bound, inter1, 'stable');
  bound2 = setdiff(P.Son2.bound, inter2, 'stable');
  
  % as the children save everything in local enumeration, we need to figure
  % out local indices
  [inter1,inter1_loc] = intersect(P.Son1.bound, inter1, 'stable');
  [inter2,inter2_loc] = intersect(P.Son2.bound, inter2, 'stable');
  [bound1,bound1_loc] = intersect(P.Son1.bound, bound1, 'stable');
  [bound2,bound2_loc] = intersect(P.Son2.bound, bound2, 'stable');
  
  % checks that everything flagged for elimination really does get eliminated
  % deactivate for performance
  if hpreconoption('checks')
    assert(isempty(union(setdiff([bound1, bound2], P.bound), setdiff(P.bound, [bound1, bound2])) ), 'Some elements in the boundary got left out');
    assert(isempty(setdiff(inter1, P.inter)), 'Found elements in inter1 that were not originally designated for elimination');
    assert(isempty(setdiff(inter2, P.inter)), 'Found elements in inter2 that were not originally designated for elimination');
    assert(isempty(setdiff(P.inter, union(inter1, inter2))), 'Found elements in P.inter that are missing in inter1 and inter2')
    assert(all(P.Son1.bound(inter1_loc) == inter1), 'Local indices of inter1 are not aligned with the global ones')
    assert(all(P.Son2.bound(inter2_loc) == inter2), 'Local indices of inter2 are not aligned with the global ones')
    assert(all(P.Son1.bound(bound1_loc) == bound1), 'Local indices of bound1 are not aligned with the global ones')
    assert(all(P.Son2.bound(bound2_loc) == bound2), 'Local indices of bound2 are not aligned with the global ones')
  end
  
  % since the ordering might have changed w.r.t the elimination tree, we
  % save the updated reordering
  P.inter = [inter1, inter2];
  P.bound = [bound1, bound2];
  
  % in order for the children to compute the Schur complements in the
  % correct order, we store the local enumeration in the children
  P.Son1.finter = inter1_loc';
  P.Son1.fbound = bound1_loc';
  P.Son2.finter = inter2_loc';
  P.Son2.fbound = bound2_loc';
end