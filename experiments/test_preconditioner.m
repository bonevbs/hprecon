rng(0)

% set up hm-toolbox and preconditoner
%hpreconoption('compression', 'hss')
%hpreconoption('compression-tolerance', 1e-3)
%hpreconoption('solve-tolerance', 1e-12)
hssoption('compression', 'qr')
hssoption('block-size', 32)
hssoption('norm', 'fro')
hssoption('threshold', 1e-12)
hodlroption('threshold',1e-12)

% load problem
load('test.mat')

% set up preconditioner
p = hprecon(elim_tree);
p.factor(A);

%X = p.apply(eye(size(A)));
%X = p.solve(A);
%X(abs(X) < 1e-3) = 0;
%perm = [p.get_interior(), p.get_boundary()];
%spy(X(perm,perm))

% set GMRES parameters
restart = 10;
tol = 1e-9;
maxit = 100 / restart;

%% run GMRES
%[x0,fl0,rr0,it0,rv0] = gmres(A,b,restart,tol,maxit);
[x1,fl1,rr1,it1,rv1] = gmres(A,b,restart,tol,maxit,@p.solve);

%fprintf('Iterations without preconditioning:              %4.0f\n', length(rv0)-1)
fprintf('Iterations with hierarchical preconditioning:    %4.0f\n', length(rv1)-1)

% plot the result
%semilogy(0:length(rv0)-1,rv0/norm(b),'-o');
%hold on
semilogy(0:length(rv1)-1,rv1/norm(p.solve(b)),'-o');
hold on
xlabel('Iteration number');
ylabel('Relative residual');
hold off