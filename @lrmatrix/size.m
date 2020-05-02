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

m = size(A.U,1);
n = size(A.V,1);

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

