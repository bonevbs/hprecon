% control random seed
rng(0)

% set GMRES parameters
restart = 10;
tol = 1e-9;
maxit = 100 / restart;

pvals = [3,4,5,6,7,8,9,10,11,12; 2,3,4,5,6,7,8,NaN,NaN,NaN; 1,2,3,4,5,6,NaN,NaN,NaN,NaN];
hvals = [16; 32; 64]*ones(1,size(pvals,2));
k = 19.5;
iters = zeros(size(pvals));
ranks = zeros(size(pvals));
dofs = zeros(size(pvals));

for ih=1:size(hvals,1)
  for ip=1:size(pvals,2)
    h = hvals(ih,ip);
    p = pvals(ih,ip);
    if isnan(p)
      iters(ih, ip) = NaN;
      ranks(ih, ip) = NaN;
      dofs(ih, ip) = NaN;
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
    
    fprintf('GMRES results with h=1/%d, p=%d and k=%4.1f\n', h, p, k)
    fprintf('Iterations with hierarchical preconditioning:    %4.0f\n', length(rv1)-1)
    
    iters(ih, ip) = length(rv1)-1;
    ranks(ih, ip) = maxrank(precon, 'complements');
    dofs(ih, ip) = size(A,2);
  end
end

%% plot iterations
figure
for ih=1:size(hvals,1)
  plot(dofs(ih,:), iters(ih,:), '-o', 'DisplayName', sprintf('h=1/%d', hvals(ih,1)))
  hold on
end
legend('Location','northeast','NumColumns',1)
%title('GMRES iterations for different wave numbers')
xlabel('# DOFs')
ylabel('# GMRES iterations')
hold off

%% plot ranks
figure
for ih=1:size(hvals,1)
  plot(dofs(ih,:), ranks(ih,:), '-o', 'DisplayName', sprintf('h=1/%d', hvals(ih,1)))
  hold on
end
legend('Location','northeast','NumColumns',1)
%title('GMRES iterations for different wave numbers')
xlabel('# DOFs')
ylabel('maximum HSS rank')
hold off

%% Output a table