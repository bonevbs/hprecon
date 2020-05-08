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

if isempty(precon_lrcompression)
	precon_lrcompression = 1;
end

% if isempty(precon_stol)
%   precon_stol = 1e-12;
% end

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
%     case 'solve-tolerance'
%       opt = precon_stol;
    case 'merging-algorithm'
      opt = precon_mergemode;
		case 'lrcompression'
			opt = precon_lrcompression;
    case 'checks'
			opt = precon_checks;
    case 'levels'
      opt = precon_levels;
    otherwise
      error('Unsupported option specified');
	end
else
	switch key
%     case 'solve-tolerance'
%       if value < 0
%         error('solve-tolerance has to be positive');
%       else
%         precon_itol = value;
%       end
    case 'merging-algorithm'
			if strcmp(value, 'direct')
        precon_mergemode = value;
      elseif strcmp(value, 'martinsson') || strcmp(value, 'direct')
        precon_mergemode = value;
        if precon_lrcompression == 0
          precon_lrcompression = 1;
          warning('lrcompression was set to 0. Activating low-rank compression of Gauss transforms.')
        end
			else
				error('Unknown merging algorithm specified.');
      end
		case 'lrcompression'
			if value == 0 || value == 1
        precon_lrcompression = value;
        if (strcmp(value, 'martinsson') || strcmp(value, 'direct')) && ~precon_lrcompression
          warning('Switching to direct merging of Schur complements.')
        end
			else
				error('Binary value required for lrcompression');
      end
    case 'checks'
      if value ~= true && value ~= false
				error('Logical value required for checks');
			else
				precon_checks = value;
      end
    case 'levels'
      if (floor(value) ~= ceil(value)) || value == 0
        error('Non-zero integer value required for number of levels');
      else
        precon_levels = value;
      end
    otherwise
      error('Unsupported option specified');
  end
end