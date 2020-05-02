function r = hprecon_solve(P,r)
  % apply recursively in post-ordering the left inverse
  r = Lsolve_rec(P, r, 1, hpreconoption('levels'));
  
  % apply diagonal blocks and the Schur complement of the top block
  r = Dsolve_rec(P, r, 1, hpreconoption('levels'));
  r(P.bound,:) = P.S \ r(P.bound,:);
  
  % apply recursively in pre-ordering the left inverse
  r = Rsolve_rec(P, r, 1, hpreconoption('levels'));
end

function r = Lsolve_rec(P, r, l, lmax)
  % process leaf nodes first
  if (l < lmax || lmax == -1) && (~isempty(P.Son1) && ~isempty(P.Son2))
    r = Lsolve_rec(P.Son1, r, l+1, lmax);
    r = Lsolve_rec(P.Son2, r, l+1, lmax);
    inter = P.inter;
  else
    inter = P.get_interior();
  end
  r(P.bound,:) = r(P.bound,:) - P.LU*(P.LV'*r(inter,:));
end

function r = Dsolve_rec(P, r, l, lmax)
	% the order doesn't matter her 
  if (l < lmax || lmax == -1) && (~isempty(P.Son1) && ~isempty(P.Son2))
    r = Dsolve_rec(P.Son1, r, l+1, lmax);
    r = Dsolve_rec(P.Son2, r, l+1, lmax);
    inter = P.inter;
  else
    inter = P.get_interior();
  end
  if ~isempty(P.Aii) && ~isstruct(P.Aii)
    r(inter,:) = P.Aii \ r(inter,:);
  elseif ~isempty(P.invAii)
    r(inter,:) = P.invAii * r(inter,:);
  elseif isstruct(P.Aii)
    inter1 = P.inter(1:size(P.Aii.A11,2));
    inter2 = P.inter(size(P.Aii.A11,2) + (1:size(P.Aii.A22,2)));
    [r1, r2] = blksolve_right(P.Aii, r(inter1,:), r(inter2,:));
    r(inter,:) = [r1; r2];
  end
end

function r = Rsolve_rec(P, r, l, lmax)
  if (l < lmax || lmax == -1) && (~isempty(P.Son1) && ~isempty(P.Son2))
    inter = P.inter;
  else
    inter = P.get_interior();
  end
  % process top nodes first
  r(inter,:) = r(inter,:) - P.RU*(P.RV'*r(P.bound,:));
  if (l < lmax || lmax == -1) && (~isempty(P.Son1) && ~isempty(P.Son2))
    r = Rsolve_rec(P.Son2, r, l+1, lmax);
    r = Rsolve_rec(P.Son1, r, l+1, lmax);
  end
end
