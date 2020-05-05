classdef lrmatrix
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
      elseif nargin == 1
        A = varargin{1};
        [U,S,V] = svd_from_random_sampling(@(x) A*x, @(x) A'*x, size(A,2), 1e-9);
        obj = lrmatrix(U*S, V);
      elseif nargin == 2
        obj.U = varargin{1};
        obj.V = varargin{2};
        assert(size(obj.U,2) == size(obj.V,2), 'Dimension mismatch.');
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