function opt = hpreconoption(key, value)
%HODLROPTION Set or get an option for nylea preconditioner.
%
% Valid options are:
%   'compression':           Possible values are none, hodlr, hss
%   'compression-tolerance': Value used for off-diagonal truncation.
%   'inversion-tolerance':   Value used for inversion of hierarchical matrices.

global precon_compression
global precon_ctol
global precon_stol
global precon_checks

if isempty(precon_compression)
	precon_compression = 'none';
end

if isempty(precon_ctol)
  precon_ctol = 1e-3;
end

if isempty(precon_stol)
  precon_stol = 1e-12;
end

if isempty(precon_checks)
  precon_checks = false;
end

if ~exist('key', 'var')
  error('Please specify a key');
end

if ~exist('value', 'var')
  switch key
    case 'compression-tolerance'
      opt = precon_ctol;
    case 'solve-tolerance'
      opt = precon_stol;
		case 'compression'
			opt = precon_compression;
    case 'checks'
			opt = precon_checks;
    otherwise
      error('Unsupported option specified');
	end
else
	switch key
    case 'compression-tolerance'
      if value < 0
        error('compression-tolerance has to be positive');
      else
        precon_ctol = value;
      end
    case 'solve-tolerance'
      if value < 0
        error('solve-tolerance has to be positive');
      else
        precon_itol = value;
      end
		case 'compression'
			if ~strcmp(value, 'none') && ~strcmp(value, 'hodlr') && ~strcmp(value, 'hss')
				error('Invalid value for compression');
			else
				precon_compression = value;
      end
    case 'checks'
      if value ~= true && value ~= false
				error('Logical value required for checks');
			else
				precon_checks = value;
      end
    otherwise
      error('Unsupported option specified');
  end
end