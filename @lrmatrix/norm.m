function nrm = norm(A,nrmtype)
if size(A.U,2) ~= size(A.V,2)
  error('uncompatible dimensions');
end
[~, RU] = qr(full(A.U),0);
[~, RV] = qr(full(A.V),0);

if ~exist('nrmtype', 'var') || ischar(nrmtype) && strcmp(nrmtype, 'fro')
  nrm = norm(RU * RV', 'fro');
else
  nrm = norm(RU * RV', nrmtype);
end

end