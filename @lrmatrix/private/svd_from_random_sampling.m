function [U, S, V] = svd_from_random_sampling(Afun, Afunt, sz, tol, reltol)

if reltol
  nrm = svds(@(x,tflag) normest_afun(Afun, Afunt, x, tflag), ...
    sz, 1, 'largest', ...
    struct('tol', 1e-2)); nrm = nrm(1);
else
  nrm = 1;
end

Q = find_range(Afun, sz(2), tol, nrm);
B = (Afunt(Q))';
[U,S,V] = svd(B,'econ');
U = Q * U;

end

function [Q, k] = find_range(Afun, n, tol, nrm)
r = 10;
Y = Afun(randn(n, r));

k = 0;
Q = [];
while norm(Y(:,k+1:end),2) > nrm * tol/(10*sqrt(2/pi))
  k = k+1;
  if k>1
    Y(:,k) = Y(:,k) - Q*Q'*Y(:,k);
  end
  q = Y(:,k)/norm(Y(:,k));
  Q = [Q, q];
  y = Afun(randn(n, 1));
  y = y - Q*Q'*y;
  Y(:,k+1:end) = Y(:,k+1:end) - q*q'*Y(:,k+1:end);
  Y = [Y,y];
end
end

function y = normest_afun(Afun, Afunt, x, tflag)
if strcmp(tflag, 'transp')
    y = Afunt(x);
else
    y = Afun(x);
end
end
