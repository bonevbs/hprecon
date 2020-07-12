classdef blkmatrix < handle
  properties
    A11
    A12
    A21
    A22
    S11
    S22
  end
  methods
    function obj = blkmatrix(varargin)
      if nargin == 0
        return;
      elseif nargin == 3
        A = varargin{1};
        mp = varargin{2};
        np = varargin{3};
        obj = blkmatrix(A(1:mp,1:np), A(1:mp,np+1:end), A(mp+1:end,1:np), A(mp+1:end,np+1:end));
      elseif nargin == 4
        obj.A11 = varargin{1};
        obj.A12 = varargin{2};
        obj.A21 = varargin{3};
        obj.A22 = varargin{4};
        if size(obj.A11,1) ~= size(obj.A12,1) || ...
            size(obj.A11,2) ~= size(obj.A21,2) || ...
            size(obj.A22,1) ~= size(obj.A21,1) || ...
            size(obj.A22,2) ~= size(obj.A12,2)
          error('Dimension mismatch. Make sure dimensions are compatible.');
        end
      end
    end
  end
end