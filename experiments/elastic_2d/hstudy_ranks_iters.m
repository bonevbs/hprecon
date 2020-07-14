% control random seed
rng(0)

% set GMRES parameters
restart = 10;
tol = 1e-9;
maxit = 100 / restart;

hvals = [16, 32, 64, 128, 256; 16, 32, 64, 128, 256; 16, 32, 64, 128, NaN];
pvals = [1; 2; 3]*ones(1,size(hvals,2));
iters = zeros(size(pvals));
ranks = zeros(size(pvals));
dofs = zeros(size(pvals));

for ih=1:size(hvals,2)
  for ip=1:size(pvals,1)
    h = hvals(ip,ih);
    p = pvals(ip,ih);
    if isnan(h)
      iters(ip, ih) = NaN;
      ranks(ip, ih) = NaN;
      dofs(ip, ih) = NaN;
      continue
    end
    
    GenMatrixElasticity2D('test.mat', h, p, 'forced 2', 10);
    load('test.mat')
    
    % set compression options
    bsz = 2*10*(p+1)*(p+2)/2;
    hpreconoption('merging-algorithm', 'martinsson')
    hpreconoption('levels', -4)
    %hpreconoption('compression-tolerance', 1e-3)
    %hpreconoption('solve-tolerance', 1e-12)
    hssoption('compression', 'svd')
    hssoption('block-size', bsz)
    hssoption('norm', 2)
    % see whether 1e-3 with high inversion accuracy works
    hssoption('threshold', 1e-6)
    % set up preconditioner
    precon = hprecon(elim_tree);
    precon.factor(A);
    
    % run GMRES
    [x1,fl1,rr1,it1,rv1] = gmres(A,b,restart,tol,maxit,@precon.solve);
    
    fprintf('GMRES results with h=1/%d, p=%d\n', h, p)
    fprintf('Iterations with hierarchical preconditioning:    %4.0f\n', length(rv1)-1)
    
    iters(ip, ih) = length(rv1)-1;
    ranks(ip, ih) = maxrank(precon, 'complements');
    dofs(ip, ih) = size(A,2);
  end
end

%% plot iterations
figure
for ip=1:size(pvals,1)
  plot(dofs(ip,:), iters(ip,:), '-o', 'DisplayName', sprintf('p=%d', pvals(ip,1)))
  hold on
end
legend('Location','northeast','NumColumns',1)
%title('GMRES iterations for different wave numbers')
xlabel('# DOFs')
ylabel('# GMRES iterations')
hold off

%% plot ranks
figure
for ip=1:size(pvals,1)
  plot(dofs(ip,:), ranks(ip,:), '-o', 'DisplayName', sprintf('p=%d', pvals(ip,1)))
  hold on
end
legend('Location','northeast','NumColumns',1)
%title('GMRES iterations for different wave numbers')
xlabel('# DOFs')
ylabel('maximum HSS rank')
hold off

%% Output a table