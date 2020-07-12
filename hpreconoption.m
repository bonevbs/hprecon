function opt = hpreconoption(key, value)
%HPRECONOPTION: Determine the 
%
% Valid options are:
%   'compression':           Possible values are none, hodlr, hss
%   'compression-tolerance': Value used for off-diagonal truncation.
%   'inversion-tolerance':   Value used for inversion of hierarchical matrices.

global precon_lrcompression
global precon_mergemode
%global precon_stol
global precon_checks
global precon_levels

if isempty(precon_mergemode)
	precon_mergemode = 'direct';
end

if isempty(precon_checks)
  precon_checks = false;
end

if isempty(precon_levels)
  precon_levels = -1;
end

if ~exist('key', 'var')
  error('Please specify a key');
end

if ~exist('value', 'var')
  switch key
    case 'merging-algorithm'
      opt = precon_mergemode;
    case 'checks'
			opt = precon_checks;
    case 'levels'
      opt = precon_levels;
    otherwise
      error('Unsupported option specified');
	end
else
	switch key
    case 'merging-algorithm'
			if strcmp(value, 'martinsson') || strcmp(value, 'direct')
        precon_mergemode = value;
			else
				error('Unknown merging algorithm specified.');
      end
    case 'checks'
      if value ~= true && value ~= false
				error('Logical value required for checks');
			else
				precon_checks = value;
      end
    case 'levels'
      if (floor(value) ~= ceil(value))
        error('Integer value required for number of levels');
      else
        precon_levels = value;
      end
    otherwise
      error('Unsupported option specified');
  end
end