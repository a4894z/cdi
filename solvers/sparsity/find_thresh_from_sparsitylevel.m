function [ thresh ] = find_thresh_from_sparsitylevel( hurr_durr, sparsity_level )

%abs_hurr_durr = abs( hurr_durr );

% sanity checking:
if ~isinteger( sparsity_level ), sparsity_level = floor( sparsity_level ); end

% if kappa = 1 (don't zero out anything):
if (size(hurr_durr,1) * size(hurr_durr,2)) <= sparsity_level, thresh = 0; return; end

% create stacked column vector:
hurr_durr_vec = hurr_durr(:);

% if kappa = 0 (i.e. kill off all pixels):
if sparsity_level <= 0, thresh = max( hurr_durr_vec ); return; end

% sort in descending order:
[ ive_been_sorted, ~ ] = sort( hurr_durr_vec, 'descend' );

% find the value to threshold corresponding to a sparsity level of kappa
thresh = ive_been_sorted( sparsity_level );

