function C = horzcat(varargin)
  if any(cellfun(@(x) size(x, 1), varargin) ~= size(varargin{1}, 1))
    error('Dimensions of arrays being concatenated are not consistent.')
  end
  C = lrmatrix();
  C.U = cell2mat(cellfun(@(x) subsref(x, struct('type','.','subs','U')), varargin, 'UniformOutput', false));
  tmp = cellfun(@(x) subsref(x, struct('type','.','subs','V')), varargin, 'UniformOutput', false);
  C.V = blkdiag(tmp{:});
end