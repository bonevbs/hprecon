S11 = hss(randn(256,256));
S22 = hss(randn(128,128));
A12 = sprand(256,128, 0.2);
A21 = sprand(128,256, 0.2);
b1 = randn(256,1);
b2 = randn(128,1);
A = [full(S11), A12; A21, full(S22)];
[z1,z2] = blocksolve_right(S11, A12, A21, S22, b1, b2);


norm(A\[b1;b2] - [z1;z2])/norm(A\[b1;b2])