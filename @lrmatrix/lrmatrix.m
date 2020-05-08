classdef lrmatrix
  % Low-rank matrices
  %
  % L = LRMATRIX(A) constructs a low-rank matrix from A using random
  %     sampling.
  %
  % L = LRMATRIX(U,V) constructs a low-rank matrix with generators U and V.
  %
  % L = LRMATRIX(Afun, Afunt, n) constructs a low-rank matrix with n
  %     columns from the handles Afun and Afunt via random sampling
  
  properties
    U
    V
  end
  
  methods
    
    function obj = lrmatrix(varargin)
      if nargin == 0
        obj.U = [];
        obj.V = [];
        return;
      elseif nargin == 1 && isa(varargin{1}, 'double')
        A = varargin{1};
        [U,S,V] = svd_from_random_sampling(@(x) A*x, @(x) A'*x, size(A), 1e-6, 0);
        obj = lrmatrix(U*S, V);
      elseif nargin == 2 && isa(varargin{1}, 'double') && isa(varargin{2}, 'double')
        obj.U = varargin{1};
        obj.V = varargin{2};
        assert(size(obj.U,2) == size(obj.V,2), 'Dimension mismatch.');
      elseif nargin == 3 && isa(varargin{1}, 'function_handle') && isa(varargin{2}, 'function_handle')
        [U,S,V] = svd_from_random_sampling(varargin{1}, varargin{2}, varargin{3}, 1e-6, 0);
        obj = lrmatrix(U*S, V);
      else
        error('Unsupported input.')
      end
    end
    
    function ind = end(obj,k,n)
      if n == 1
        ind = prod(size(obj));
      else
        ind = size(obj,k);
      end
    end
    
    function rk = rank(obj)
      rk = size(obj.U,2);
    end
    
    function obj = ctranspose(A)
      obj = lrmatrix(A.V, A.U);
    end
    
    function A = full(obj)
      A = obj.U*obj.V';
    end
  end
end