function [ sparsity_level ] = find_sparsitylevel_from_thresh( hurr_durr, thresh, N )

hurr_durr = ( hurr_durr > thresh );

sparsity_level = sum(sum( hurr_durr )) / N.NcNr;