function C = vertcat(varargin)
  if any(cellfun(@(x) size(x, 2), varargin) ~= size(varargin{1}, 2))
    error('Dimensions of arrays being concatenated are not consistent.')
  end
  C = lrmatrix();
  C.V = cell2mat(cellfun(@(x) subsref(x, struct('type','.','subs','V')), varargin, 'UniformOutput', false));
  tmp = cellfun(@(x) subsref(x, struct('type','.','subs','U')), varargin, 'UniformOutput', false);
  C.U = blkdiag(tmp{:});
end