function B = hss_martinsson_adaptive(obj, Afun, Afunt, Aeval, m, n, kest)
%HSS_MARTINSSON_ADAPTIVE Build the HSS representation of a matrix
% 			 using mat-vec multiplication with random (block) vectors
%			   and access to (block) diagonal entries [1].
%
% [1] Martinsson, Per-Gunnar. "A fast randomized algorithm for computing a
%     hierarchically semiseparable representation of a matrix." SIAM
%     Journal on Matrix Analysis and Applications 32.4 (2011): 1251-1274

reltol = false;

tol = hssoption('threshold');
nrm_type = hssoption('norm');

B = obj;

if B.topnode == 1 && B.leafnode == 1
  [m, n] = size(B);
  B.D = Aeval((1:m).', 1:n);
  return
end

failed = true;

if ~exist('kest', 'var')
  k = 10;
else
  k = kest;
end

a = 10;
%bs = ceil(0.01*n);
bs = 20;

nrm = svds(@(x,t) normest_afun(Afun, Afunt, x, t), [m n], 1, 'largest', struct('tol', 1e-2));
nrm = nrm(1);

Ocol = randn(n, k + a);
Scol = Afun(Ocol);
Orow = randn(m, k + a);
Srow = Afunt(Orow);

B = hss_extract_diagonal_rec(B, Aeval, 0, 0);

while failed && k < n
  %bs = k;
  
  Ocol_test = randn(n, bs);
  Scol_test = Afun(Ocol_test);
  Orow_test = randn(m, bs);
  Srow_test = Afunt(Orow_test);
  [B, ~, ~, ~, ~, ~, ~, ~, ~] = hss_martinsson_rec(B, Aeval, Scol, Srow, ...
    Ocol, Orow, 0, 0, k);
  
  if reltol
    nrm_est = sqrt(1/bs)*norm(Scol_test - B*Ocol_test, 'fro')/nrm;
  else
    nrm_est = sqrt(1/bs)*norm(Scol_test - B*Ocol_test, 'fro');
  end
  failed = (nrm_est > tol);
  
  if failed
    %fprintf('HSS_MARTINSSON_ADAPTIVE :: Enlarging sampling space to %d\n', k+bs);
    Ocol = [ Ocol, Ocol_test ];
    Scol = [ Scol, Scol_test ];
    Orow = [ Orow, Orow_test ];
    Srow = [ Srow, Srow_test ];
    k = k+bs;
  end
end

% Remove any temporary data that we might have stored in the leafnodes B12
% anb B21 -- these are not part of the final HSS data structure.
B = clean_structure(B);

if failed
  error('Failed to converge.');
end

% [ hssrank(B), k ]
% B = hss_proper(B);
%B = compress(B, tol);
% hssrank(B)

end

function [B, Scol, Srow, Ocol, Orow, Jcol, Jrow, U, V] =...
  hss_martinsson_rec(B, Aeval, Scol, Srow, Ocol, Orow, row, col, k)

if B.leafnode == 1
  [m, n] = size(B.D);
  
  % substract diagonal block
  %B.D = Aeval(row + 1: row + m, col + 1:col + n);
  Scol = Scol - B.D * Ocol;
  Srow = Srow - B.D' * Orow;
  
  % take care of the columns
  %Qcol = colspan(Scol, k);
  %[Xcol, Jcol] = interpolate(Qcol');
  [Xcol, Jcol] = interpolate(Scol');
  B.U = Xcol';
  Scol = Scol(Jcol, :);
  U = Xcol';
  B.B12 = Jcol;
  Jcol = row + Jcol;
  
  % take care of the rows
  %Qrow = colspan(Srow, k);
  %[Xrow, Jrow] = interpolate(Qrow');
  [Xrow, Jrow] = interpolate(Srow');
  B.V = Xrow';
  Srow = Srow(Jrow, :);
  V = Xrow';
  B.B21 = Jrow;
  Jrow = col + Jrow;
else
  [B.A11, Scol1, Srow1, Ocol1, Orow1, Jcol1, Jrow1, U1, V1]  = hss_martinsson_rec(B.A11, Aeval, Scol(1:B.ml, :), Srow(1:B.nl, :), Ocol(1:B.nl, :), Orow(1:B.ml, :), row, col, k);
  [B.A22, Scol2, Srow2, Ocol2, Orow2, Jcol2, Jrow2, U2, V2]  = hss_martinsson_rec(B.A22, Aeval, Scol(B.ml + 1:end, :), Srow(B.nl + 1:end,:), Ocol(B.nl + 1:end, :), Orow(B.ml + 1:end, :), row + B.ml, col + B.nl, k);
  
  % update the sampling matrix based on the extracted generators
  Ocol2 = V2' * Ocol2;
  Ocol1 = V1' * Ocol1;
  Orow2 = U2' * Orow2;
  Orow1 = U1' * Orow1;
  
  % step 1 in the algorithm: look only at extracted rows/cols
  Jcol = [Jcol1, Jcol2]; Jrow = [Jrow1, Jrow2];
  Ocol = [Ocol1; Ocol2]; Orow = [Orow1; Orow2];
  
  % extract the correct generator blocks
  B.B12 = Aeval(Jcol1, Jrow2);
  B.B21 = Aeval(Jcol2, Jrow1);
  
  % subtract the diagonal block
  Scol = [Scol1 - B.B12  * Ocol2;  Scol2 - B.B21  * Ocol1 ];
  Srow = [Srow1 - B.B21' * Orow2;  Srow2 - B.B12' * Orow1 ];
  
  % if this is the top node we can stop here
  if B.topnode
    Scol = []; Srow = []; Ocol = []; Orow = [];
    Jcol = []; Jrow = []; U = []; V = [];
    failed = false;
  else
    % take care of the columns
    %Qcol = colspan(Scol, k);
    %[Xcol, Jcolloc] = interpolate(Qcol');
    [Xcol, Jcolloc] = interpolate(Scol');
    B.Rl = Xcol(:, 1:size(Scol1, 1))';
    B.Rr = Xcol(:, size(Scol1, 1)+1:end)';
    Scol = Scol(Jcolloc, :);
    Jcol = Jcol(Jcolloc);
    U = [ B.Rl ; B.Rr ];
    
    % take care of the rows
    %Qrow = colspan(Srow, k);
    %[Xrow, Jrowloc] = interpolate(Qrow');
    [Xrow, Jrowloc] = interpolate(Srow');
    B.Wl = Xrow(:, 1:size(Srow1, 1))';
    B.Wr = Xrow(:, size(Srow1, 1)+1:end)';
    Srow = Srow(Jrowloc, :);
    Jrow = Jrow(Jrowloc);
    V = [ B.Wl ; B.Wr ];
    
    B.U = Jcolloc;
    B.V = Jrowloc;
  end
end

end

function B =...
  hss_extract_diagonal_rec(B, Aeval, row, col)

if B.leafnode == 1
  [m, n] = size(B.D);
  % extract the diagonal block
  B.D = Aeval(row + 1: row + m, col + 1:col + n);  
else
  B.A11 = hss_extract_diagonal_rec(B.A11, Aeval, row, col);
  B.A22 = hss_extract_diagonal_rec(B.A22, Aeval, row + B.ml, col + B.nl);
end

end

function Q = colspan(S, k)
use_qr = false;

if use_qr
  [Q, R, ~] = qr(S, 0);
  k = min(k, size(Q,2))
  Q = Q(:,1:k);
else
  [Q, ~, ~] = svd(S);
  k = min(k, size(Q,2));
  Q = Q(:,1:k);
end
end

% axuiliary routine to clean up HSS structure, should become obsolete
function B = clean_structure(B)
%CLEAN_STRUCTURE Clean the HSS structure from temporary data.

if B.leafnode
  B.B12 = [];
  B.B21 = [];
else
  B.A11 = clean_structure(B.A11);
  B.A22 = clean_structure(B.A22);
end

end

% interpolative decomposition
function [X, J] = interpolate(A, tol, reltol)

if isempty(A)
  J = []; X = A;
  return;
end
if ~exist('reltol', 'var')
  reltol = false;
end
if ~exist('tol', 'var')
  tol = eps;
end

[~, R, p] = qr(A, 'vector');
for j = 1:length(p) % invert permutation
  ip(p(j))= j;
end
if reltol
  rk = sum(abs(diag(R)) > tol);
else
  rk = sum(abs(diag(R)) > tol * abs(R(1,1)));
end
J = p(1:rk);
X = R(1:rk, 1:rk)\ R(1:rk,ip);
end

% auxiliary function for norm estimation
function y = normest_afun(Afun, Afunt, x, transp)

if strcmp(transp, 'transp')
  y = Afunt(x);
else
  y = Afun(x);
end

end

