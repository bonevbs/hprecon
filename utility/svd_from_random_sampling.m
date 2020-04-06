function [U, S, V] = svd_from_random_sampling(Afun, Afunt, n, tol)
  k = 10;
  bs = 5;
  
  O1 = randn(n, k);
  S1 = Afun(O1);
  
  success = false;
  i = 1; maxit = 100;
	while ~success && i < maxit
    [Q,~] = qr(S1,0);
    O2 = randn(n, bs);
    S2 = Afun(O2);

    %norm(S2 - Q*Q'*S2,'fro')/norm(S2,'fro');
    success = norm(S2 - Q*Q'*S2,'fro') < tol * norm(S2,'fro');
    
    if ~success
      O1 = [O1, O2];
      S1 = [S1, S2];
      k = k + bs;
    end
    i = i+1;
  end
  
  assert(success, 'Approximate SVD did not converge!');
  
  B = (Afunt(Q))';
  [U,S,V] = svd(B,'econ');
  U = Q * U;
  
end