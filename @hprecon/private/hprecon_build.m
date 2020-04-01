function P = hprecon_build(elim_tree)
  root = find(elim_tree.fathers == -1);
  assert(isscalar(root), 'More than one root has been found.')
  [P,~] = hprecon_build_rec(elim_tree, root);
end

function [P, ncreated] = hprecon_build_rec(elim_tree, i)
  P = hprecon();
  n1 = 0; n2 =0;
  
  if (elim_tree.lsons(i) ~= -1)
    [P.Son1, n1] = hprecon_build_rec(elim_tree, elim_tree.lsons(i));
  end
  if (elim_tree.rsons(i) ~= -1)
    [P.Son2, n2] = hprecon_build_rec(elim_tree, elim_tree.rsons(i));
  end
  
  P.desc = n1 + n2;
  ncreated = P.desc + 1;
  ninter = elim_tree.ninter(i);
  nbound = elim_tree.nbound(i);
  P.inter = elim_tree.inter(1:ninter,i)';
  P.bound = elim_tree.bound(1:nbound,i)';
end