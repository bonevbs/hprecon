function maxl = maxlevel(obj, varargin)
% MAXRANK(p): traverses the preconditioner tree structure to find the
% number of levels in the tree structure
maxl = maxlevel_rec(obj, 1, 1);
end

function maxl = maxlevel_rec(obj, l, maxl)
  if (isempty(obj.Son1) && isempty(obj.Son2))
    maxl = max(l, maxl);
  elseif ~isempty(obj.Son1) && ~isempty(obj.Son2)
    maxl = maxlevel_rec(obj.Son1, l+1, maxl);
    maxl = maxlevel_rec(obj.Son2, l+1, maxl);
  else
    error('MAXLEVEL: Something went wrong. This is not supported.')
  end
end