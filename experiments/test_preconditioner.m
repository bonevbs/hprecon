rng(0)

% set up hm-toolbox and preconditoner
preconoption('compression', 'hss')
preconoption('compression-tolerance', 1e-3)
preconoption('solve-tolerance', 1e-12)
hssoption('compression', 'svd')
hssoption('block-size', 32)
hssoption('norm', 'fro')
hssoption('threshold', 1e-3)

% load problem
load('test.mat')

% set up preconditioner
p = hprecon(elim_tree);
p.factor(A);
