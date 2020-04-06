rng(0);

k = 50;
U = randn(300,k);
V = randn(200,k);
[U, ~] = qr(U,0);
[V, ~] = qr(V,0);

A = U*(diag(exp(-1*(1:k)))*V');

[L,S,R] = svd_from_random_sampling(@(x) A*x, @(x) A'*x, size(A,2), 1e-1);
fprintf("Detected a rank of %d\n", size(S,2));