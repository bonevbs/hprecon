classdef blkmatrix
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
        obj.A11 = [];
        obj.A12 = [];
        obj.A21 = [];
        obj.A22 = [];
        obj.S11 = [];
        obj.S22 = [];
        return;
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