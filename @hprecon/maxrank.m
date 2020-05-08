function maxk = maxrank(p, varargin)
% MAXRANK(p): traverses the preconditioner tree structure to find the
% maximum rank
relrank = 0;
mode = 'both';
if nargin > 1
  if strcmp(varargin{1}, 'complements') || strcmp(varargin{1}, 'transforms') || strcmp(varargin{1}, 'both')
    mode = varargin{1};
  else
    error('Unknown mode.')
  end
end
if nargin > 2
  if strcmp(varargin{2}, 'relative')
    relrank = 1;
  elseif strcmp(varargin{2}, 'absolute')
    relrank = 0;
  else
    error('Specify whether ranks should be computed in a relative way.')
  end
end


maxk =  maxrank_rec(p, 0, mode, relrank, 1);

end

function [maxk, i] = maxrank_rec(p, maxk, mode, relrank, i)

if ~isempty(p.Son1)
  [maxk, i] = maxrank_rec(p.Son1, maxk, mode, relrank, i);
end
if ~isempty(p.Son2)
  [maxk, i] = maxrank_rec(p.Son2, maxk, mode, relrank, i);
end

if ~isempty(p.S)
  if strcmp(mode, 'complements') || strcmp(mode, 'both') 
    if ~isa(p.S, 'hss')
      warning('Schur complement is not HSS. Ignoring the hss-ranks.')
    else
      k = hssrank(p.S);
      if relrank
        k = k/size(p.S,2);
      end
      maxk = max(k, maxk);
      fprintf('Rank of Schur complement #%d: %d \n', i, k);
    end
  end
  if strcmp(mode, 'transforms') || strcmp(mode, 'both')
    if isa(p.L, 'lrmatrix')
      if relrank
        lk = rank(p.L)/min(size(p.L));
      else
        lk = rank(p.L);
      end
      maxk = max(lk, maxk);
      fprintf('Rank of left Gauss transform #%d: %d \n', i, lk);
    end
    if isa(p.R, 'lrmatrix')
      if relrank
        rk = rank(p.R)/min(size(p.R));
      else
        rk = rank(p.R);
      end
      maxk = max(rk, maxk);
      fprintf('Rank of right Gauss transform #%d: %d \n', i, rk);
    end
  end
  i = i + 1;
end

end