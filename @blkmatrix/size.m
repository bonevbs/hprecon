function [sz, sz2] = size(A, idx)
%SIZE Size of an HSS matrix.

if nargout == 2
  if exist('idx', 'var')
    error('The dimension cannot be specified if two outputs are expected');
  end
  
  sz  = size(A, 1);
  sz2 = size(A, 2);
  
  return;
end

if size(obj.A11,1) ~= size(obj.A12,1) || ...
    size(obj.A11,2) ~= size(obj.A21,2) || ...
    size(obj.A22,1) ~= size(obj.A21,1) || ...
    size(obj.A22,2) ~= size(obj.A12,2)
  error('Dimension mismatch. Something went wrong');
end

m = size(A.A11,1) + size(A.A21,1);
n = size(A.A11,2) + size(A.A11,2);

if exist('idx', 'var')
  switch idx
    case 1
      sz = m;
    case 2
      sz = n;
  end
else
  sz = [m n];
end

end

