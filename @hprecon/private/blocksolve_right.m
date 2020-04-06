function [x1, x2] = blocksolve_right(A11, A12, A21, A22, b1, b2)
  % A11 and A22 must be in hss form
  x1 = A11\b1;
  x2 = b2 - A21*x1;
  S22 = form_hss_complement(A11, A12, A21, A22);
  %fprintf("BLOCKSOLVE_RIGHT: rel. error in the S22 Schur complement: %e\n",...
  %  norm(A22 - A21*(A11\A12) - S22, 'fro')/norm(A22 - A21*(A11\A12), 'fro'));
  x2 = S22\x2;
  x1 = x1 - A11\(A12*x2);
end

function S22 = form_hss_complement(A11, A12, A21, A22)
  % use nested functions in order to adjust the tolerance for inversion
  function z = Sfun(x)
    tol = hssoption('threshold');
    hssoption('threshold', 1e-12);
    z = A22*x - A21*(A11\(A12*x));
    hssoption('threshold', tol);
  end
  
  function z = Sfunt(x)
    tol = hssoption('threshold');
    hssoption('threshold', 1e-12);
    z = A22'*x - A12'*((A11')\(A21'*x));
    hssoption('threshold', tol);
  end

  function S = Seval(i,j)
    tol = hssoption('threshold');
    hssoption('threshold', 1e-12);
    S = get(A22, i, j) - A21(i,:)*(A11\A12(:,j));
    hssoption('threshold', tol);
  end
  
  %Sfun = @(x) A22*x - A21*(A11\(A12*x));
  %Sfunt = @(x) A22'*x - A12'*((A11')\(A21'*x));
  % the evaluation step is currently quadratic in complexiy.
  % a more efficient implementation would take the relevant indices of
  % inv(A11) and only look at the submatrix of A11
  %Seval = @(i,j) get(A22, i, j) - A21(i,:)*(A11\A12(:,j));
  
  S22 = hss('handle', @(x) Sfun(x), @(x) Sfunt(x), @(i,j) Seval(i,j), size(A22,1), size(A22,2));
end