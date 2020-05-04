classdef hprecon < handle
  properties
    % tree structure
    desc
    Son1
    Son2
    % index sets corresponding to global nodes
    inter
    bound
    % index sets corresponding to local enumeration of the parent interior
    % and boundary nodes
    finter
    fbound
    % interior block and Schur complement in HSS format
    invAii
    Aii
    S
    % factor for the left Gauss inverse
    L
    % factors for the right Gauss inverse
    R
  end
  
  methods
    
    % constructor
    function obj = hprecon(elim_tree)
      if nargin == 0
        obj.desc   = 0;
        obj.Son1   = [];
        obj.Son2   = [];
        obj.inter  = [];
        obj.bound  = [];
        obj.finter = [];
        obj.fbound = [];
        return;
      elseif nargin == 1
        obj = hprecon_build(elim_tree);
      end
    end
    
    % compute factorization
    function factor(obj, A)
      if strcmp(hpreconoption('compression'),'both')
        hprecon_fact_lr(obj,A);
      elseif strcmp(hpreconoption('compression'),'hss')
        hprecon_fact_full(obj,A);
      else
        error('Currently unsupported!\n');
      end
    end
    
    
    % apply preconditioner
    function x = solve(obj, r)
      x = r;
      x = hprecon_solve(obj, x);
    end
    
    % spits out all DOFs that were eliminated
    function interior = get_interior(obj)
      interior = [];
      interior = get_interior_rec(obj, interior);
    end
    
    % spits out bounding box
    function boundary = get_boundary(obj)
      boundary = obj.bound;
    end
    
    % routine to check the complements
    function check_complements(obj, A)
      check_complements_rec(obj, A, 1);
    end
  end
end