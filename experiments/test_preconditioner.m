rng(0)

% parameters of the discretization
pdeg = 3;
hinv = 256;
bsz = 2*10*(pdeg+1)*(pdeg+2)/2;

%GenMatrixPoisson2D('test.mat', hinv,pdeg, 10)
%GenMatrixPoissonCG2D('test.mat', hinv,pdeg, 25)
%GenMatrixHighContrast2D('test.mat', 64, 2, 25, 1, 0);
%GenMatrixHelmholtz2D('test.mat', hinv, pdeg, 10, 19.5);
GenMatrixElasticity2D('test.mat', hinv, pdeg, 'forced 2', 10);
% load problem
load('test.mat')

% set up hm-toolbox and preconditoner
%hpreconoption('lrcompression', 0)
hpreconoption('merging-algorithm', 'martinsson')
hpreconoption('levels', -4)
%hpreconoption('compression-tolerance', 1e-3)
%hpreconoption('solve-tolerance', 1e-12)
hssoption('compression', 'svd')
hssoption('block-size', bsz)
hssoption('norm', 2)
% see whether 1e-3 with high inversion accuracy works
hssoption('threshold', 1e-6)
% implement option to force same depth on all bottom level Schur complement
% hierarchies

% set up preconditioner
p = hprecon(elim_tree);
p.factor(A)

%X = p.solve(eye(size(A)));
%X = p.solve(A);
%eps * cond(X)
%X(abs(X) < 1e-3) = 0;
%perm = [p.get_interior(), p.get_boundary()];
%spy(X(perm,perm))

% set GMRES parameters
restart = 10;
tol = 1e-9;
maxit = 100 / restart;

%fprintf('Condition number of matrix:                   %4.0f\n', cond(full(A)))
%fprintf('Condition number of preconditioned matrix:    %4.0f\n', cond(full(p.solve(A))))

%e = eig(full(p.solve(A)))
%figure
%plot(e,'rx');

%% run GMRES
%[x0,fl0,rr0,it0,rv0] = gmres(A,b,restart,tol,maxit);
[x1,fl1,rr1,it1,rv1] = gmres(A,b,restart,tol,maxit,@p.solve);
%[x2,fl2,rr2,it2,rv2] = gmres(A,b,restart,tol,maxit,@(r) X*r);

%fprintf('Iterations without preconditioning:              %4.0f\n', length(rv0)-1)
fprintf('Iterations with hierarchical preconditioning:    %4.0f\n', length(rv1)-1)
%fprintf('Iterations with matrix representation       :    %4.0f\n', length(rv2)-1)

% plot the result
figure
%semilogy(0:length(rv0)-1,rv0/norm(b),'-o');
%hold on
semilogy(0:length(rv1)-1,rv1/norm(p.solve(b)),'-o');
hold on
xlabel('Iteration number');
ylabel('Relative residual');
hold off