rng(0)

%GenMatrixHighContrast2D('test.mat', 64, 2, 100, 1, 0);
GenMatrixHelmholtz2D('test.mat', 64, 2, 25, 15);

% set up hm-toolbox and preconditoner
hpreconoption('compression', 'both')
hpreconoption('levels', 8)
%hpreconoption('compression-tolerance', 1e-3)
%hpreconoption('solve-tolerance', 1e-12)
hssoption('compression', 'svd')
hssoption('block-size', 30)
hssoption('norm', 2)
hssoption('threshold', 1e-6)


% load problem
load('test.mat')

% set up preconditioner
p = hprecon(elim_tree);
p.factor(A)

%X = p.solve(eye(size(A)));
%X = p.solve(A);
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
%[x1,fl1,rr1,it1,rv1] = gmres(A,b,restart,tol,maxit,@(r) X*r);

%fprintf('Iterations without preconditioning:              %4.0f\n', length(rv0)-1)
fprintf('Iterations with hierarchical preconditioning:    %4.0f\n', length(rv1)-1)

% plot the result
figure
%semilogy(0:length(rv0)-1,rv0/norm(b),'-o');
%hold on
semilogy(0:length(rv1)-1,rv1/norm(p.solve(b)),'-o');
hold on
xlabel('Iteration number');
ylabel('Relative residual');
hold off