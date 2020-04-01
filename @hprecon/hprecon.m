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
    Aii
    S
    % factor for the left Gauss inverse
    LU
    LV
    % factors for the right Gauss inverse
    RU
    RV
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
        obj = build_hprecon(elim_tree);
      end
    end
    
    % compute factorization
    function factor(obj, A)
      symbolic_factorization_rec(obj,1);
      compute_factorization_rec(obj,A,1);
    end
    
    
    % apply preconditioner
    function x = apply(obj, r)
      x = r;
      x = apply_hprecon(obj, x);
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