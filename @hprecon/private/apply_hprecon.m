function r = apply_hprecon(P,r)
  % apply recursively in post-ordering the left inverse
  r = apply_CinvA_rec(P, r);
  
  % apply diagonal blocks and the Schur complement of the top block
  r = apply_invA_rec(P, r);
  r(P.bound,:) = P.S\r(P.bound,:);
  
  % apply recursively in pre-ordering the left inverse
  r = apply_invAB_rec(P, r);
end

function r = apply_CinvA_rec(P, r)
  % process leaf nodes first
  if (~isempty(P.Son1) && ~isempty(P.Son2))
    r = apply_CinvA_rec(P.Son1, r);
    r = apply_CinvA_rec(P.Son2, r);
  end
  r(P.bound,:) = r(P.bound,:) - P.CinvA*r(P.inter,:);
end

function r = apply_invA_rec(P, r)
	% the order doesn't matter here
  if (~isempty(P.Son1) && ~isempty(P.Son2))
    r = apply_invA_rec(P.Son1, r);
    r = apply_invA_rec(P.Son2, r);
  end
  r(P.inter,:) = P.invA*r(P.inter,:);
end

function r = apply_invAB_rec(P, r)
  % process top nodes first
  r(P.inter,:) = r(P.inter,:) - P.invAB*r(P.bound,:);
  if (~isempty(P.Son1) && ~isempty(P.Son2))
    r = apply_invAB_rec(P.Son2, r);
    r = apply_invAB_rec(P.Son1, r);
  end
end
