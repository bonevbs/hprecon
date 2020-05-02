function [x1, x2] = blksolve_right(blkA, b1, b2)
  A11 = blkA.A11; A12 = blkA.A12; A21 = blkA.A21; A22 = blkA.A22;
  % A11 and A22 must be in hss form
  x1 = A11\b1;
  x2 = b2 - A21*x1;
  %hssoption('threshold', tol);
  %tol = hssoption('threshold');
  %hssoption('threshold', 1e-3*tol);
  %S22 = form_hss_complement(A11, A12, A21, A22);
  %hssoption('threshold', tol);
  S22d = hss(A22 - A21*(A11\A12));
  %S22 = (A22 - A21*(A11\A12));
  %S22 = hss('handle', @(x) S22*x, @(x) S22'*x, @(i,j) S22(i,j), size(A22,1), size(A22,2));
  %hssoption('threshold', tol);
  %fprintf("BLOCKSOLVE_RIGHT: hss rank of the S22 Schur complement: %d\n", hssrank(S22));
%   fprintf("BLOCKSOLVE_RIGHT: rel. error of the Schur complement: %e\n",...
%     norm((A22 - A21*(A11\A12)) - full(S22))/norm(A22 - A21*(A11\A12)));
%   fprintf("BLOCKSOLVE_RIGHT: rel. error of Schur complement with direct compression: %e\n",...
%     norm((A22 - A21*(A11\A12)) - full(S22d))/norm(A22 - A21*(A11\A12)));
%   fprintf("BLOCKSOLVE_RIGHT: rel. error in the inverse of the Schur complement: %e\n",...
%     norm(inv(A22 - A21*(A11\A12)) - inv(S22))/norm(inv(A22 - A21*(A11\A12))));
%   fprintf("BLOCKSOLVE_RIGHT: rel. error in the inverse of the Schur complement with direct compression: %e\n",...
%     norm(inv(A22 - A21*(A11\A12)) - inv(S22d))/norm(inv(A22 - A21*(A11\A12))));
  %fprintf("BLOCKSOLVE_RIGHT: rel. error of the solution: %e\n",...
  %  norm((A22 - A21*(A11\A12))\x2 - S22\x2)/norm((A22 - A21*(A11\A12))\x2));
  %if norm((A22 - A21*(A11\A12))\x2 - S22\x2)/norm((A22 - A21*(A11\A12))\x2) > 1e-1
  %  warning("BLOCKSOLVE_RIGHT: big error");
  %end
  %hssoption('threshold', 1e-10);
  %norm(x2)
  %norm(x2 - (A22 - A21*(A11\A12))*(S22\x2))/norm(x2)
  %norm(x2 - (A22 - A21*(A11\A12))*(S22d\x2))/norm(x2)
%   e1 = norm(inv(A22 - A21*(A11\A12)) - inv(S22));
%   e2 = norm(inv(A22 - A21*(A11\A12)) - inv(S22d));
%   if (e1 > e2 && e1 > 0.001)
%     norm((A22 - A21*(A11\A12))\x2 - S22\x2);
%     norm((A22 - A21*(A11\A12))\x2 - S22d\x2);
%     
%     % nada
%     %warning('Replacing S22')
%     %S = A22 - A21*(A11\A12);
%     %S22 = hss('handle', @(x) S*x, @(x) S'*x, @(i,j) S(i,j), size(S,1), size(S,2));
%     %S22 = S22d;
%     
%     %S22.B12 = S22d.B12;
%     %S22.B21 = S22d.B21;
%     %S22.A11.U = S22d.A11.U;
%     %S22.A11.V = S22d.A11.V;
%     %S22.A22.U = S22d.A22.U;
%     %S22.A22.V = S22d.A22.V;
%     
%     %norm(S22.A11.U*S22.B12*S22.A22.V' - S22d.A11.U*S22d.B12*S22d.A22.V')
%     tol = hssoption('threshold');
%     [U,S,V] = svd(S22.A11.U*S22.B12*S22.A22.V');
%     rk = min(find(cumsum(diag(S)'.^2) >= (1-tol^2)*sum(diag(S)'.^2)))
%     rk = sum(diag(S)' > tol * S(1,1))
%     
%     
%     %S22.A11.D = S22d.A11.D;
%     %S22.A22.D = S22d.A22.D;
%     %S22.B12 = S22.B12(1:size(S22d.B12,1),1:size(S22d.B12,2));
%     %S22.B21 = S22.B21(1:size(S22d.B21,1),1:size(S22d.B21,2));
%     %S22.A11.U = S22.A11.U(:, 1:size(S22d.A11.U,2));
%     %S22.A11.V = S22.A11.V(:, 1:size(S22d.A11.V,2));
%     %S22.A22.U = S22.A22.U(:, 1:size(S22d.A22.U,2));
%     %S22.A22.V = S22.A22.V(:, 1:size(S22d.A22.V,2));
%   end
  %x2 = (A22 - A21*(A11\A12))\x2;
  %F = ulv(S22);
  %x2 = ulv_solve(F, x2);
  %tol = hssoption('threshold');
  %hssoption('threshold', 1e-12);
  x2 = S22d\x2;
  x1 = x1 - A11\(A12*x2);
  %hssoption('threshold', tol);
end

function S22 = form_hss_complement(A11, A12, A21, A22)
  % use nested functions in order to adjust the tolerance for inversion
  function z = Sfun(x)
    z = A22*x - A21*(A11\(A12*x));
  end
  
  function z = Sfunt(x)
    z = (A22')*x - (A12')*((A11')\((A21')*x));
  end

  function S = Seval(i,j)
    %S = A22 - A21*(A11\A12);
    %S = S(i,j);
    S = get(A22, i, j) - A21(i,:)*(A11\A12(:,j));
  end
  
  %Sfun = @(x) A22*x - A21*(A11\(A12*x));
  %Sfunt = @(x) A22'*x - A12'*((A11')\(A21'*x));
  % the evaluation step is currently quadratic in complexiy.
  % a more efficient implementation would take the relevant indices of
  % inv(A11) and only look at the submatrix of A11
  %Seval = @(i,j) get(A22, i, j) - A21(i,:)*(A11\A12(:,j));
  S22 = hss('handle', @(x) Sfun(x), @(x) Sfunt(x), @(i,j) Seval(i,j), size(A22,1), size(A22,2));
end