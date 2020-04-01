function r = hprecon_solve(P,r)
  % apply recursively in post-ordering the left inverse
  r = Lsolve_rec(P, r);
  
  % apply diagonal blocks and the Schur complement of the top block
  r = Dsolve_rec(P, r);
  r(P.bound,:) = P.S \ r(P.bound,:);
  
  % apply recursively in pre-ordering the left inverse
  r = Rsolve_rec(P, r);
end

function r = Lsolve_rec(P, r)
  % process leaf nodes first
  if (~isempty(P.Son1) && ~isempty(P.Son2))
    r = Lsolve_rec(P.Son1, r);
    r = Lsolve_rec(P.Son2, r);
  end
  r(P.bound,:) = r(P.bound,:) - P.LU*(P.LV'*r(P.inter,:));
end

function r = Dsolve_rec(P, r)
	% the order doesn't matter here
  if (~isempty(P.Son1) && ~isempty(P.Son2))
    r = Dsolve_rec(P.Son1, r);
    r = Dsolve_rec(P.Son2, r);
  end
  r(P.inter,:) = P.Aii \ r(P.inter,:);
end

function r = Rsolve_rec(P, r)
  % process top nodes first
  r(P.inter,:) = r(P.inter,:) - P.RU*(P.RV'*r(P.bound,:));
  if (~isempty(P.Son1) && ~isempty(P.Son2))
    r = Rsolve_rec(P.Son2, r);
    r = Rsolve_rec(P.Son1, r);
  end
end
