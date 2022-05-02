function [ sol, expt ] = define_initial_SCPM( sol, expt )
    
if sol.init.use_prev_SCPM
    
    %====================================
    % load previous SCPMs here if desired
    %====================================

%     old_scpm = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/probes/L0105_to_L0113_combined_512x512_12Jan2022_t130245_it15001.mat', 'probe' );
%     old_scpm = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/probes/L0105_to_L0113_combined_512x512_19Apr2022_t195009_epoch2500.mat', 'probe' );
%     old_scpm = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/probes/L0274_to_L0280_combined_512x1024_20Apr2022_t022421_epoch3000.mat', 'probe' );
%     old_scpm = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/results/L0020_to_L0032_combined_512x1024_21Apr2022_t180119_epoch1000.mat', 'probe' );
%     old_scpm = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/results/L0020_to_L0032_combined_512x1024_21Apr2022_t224635_epoch3000.mat', 'probe' );
%     old_scpm = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/results/L0117_to_L0120_combined_768x1024_24Apr2022_t145039_epoch3000.mat', 'probe' );
   old_scpm = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/results/L0117_to_L0120_combined_768x1024_25Apr2022_t022803_epoch4000.mat', 'probe' );
    
    sol.probe.phi     = old_scpm.probe.phi;
    sol.probe.scpm    = old_scpm.probe.scpm;
    sol.probe.support = old_scpm.probe.support;

    %========
    
    sol.probe.phi = [];
    
%     sol.probe.phi( :, :, 5 ) = old_scpm.probe.phi( :, :, 8 );
%     sol.probe.phi( :, :, 4 ) = old_scpm.probe.phi( :, :, 7 );
    sol.probe.phi( :, :, 3 ) = old_scpm.probe.phi( :, :, 5 );
    sol.probe.phi( :, :, 2 ) = old_scpm.probe.phi( :, :, 4 ); 
    sol.probe.phi( :, :, 1 ) = old_scpm.probe.phi( :, :, 3 );
    
    sol.probe.scpm.N = size( sol.probe.phi, 3 );

    %====================================================
    % if necessary, resize the probe modes we just loaded
    %====================================================
    
    if size( sol.probe.phi, 1 ) ~= expt.sz.r || size( sol.probe.phi, 2 ) ~= expt.sz.c
        
        tmp0 = zeros( expt.sz.sz, 'single' );
        
        for pp = 1 : sol.probe.scpm.N 

            tmp0( :, :, pp ) = imresize( sol.probe.phi( :, :, pp ), expt.sz.sz );

        end
        
        sol.probe.phi = tmp0;
        
    end
    
    %=====================
    % modify current SCPMs
    %=====================    
    
%     for pp = [ 1 ]
% 
%         probe.fwhm.r = single( 0.8e-6 / sol.csys.z2.dLy );
%         probe.fwhm.c = single( 0.8e-6 / sol.csys.z2.dLx );
% 
%         probe.fwhm.r = (( 2 * rand - 1 ) * 0.05 + 1.0 ) * probe.fwhm.r;
%         probe.fwhm.c = (( 2 * rand - 1 ) * 0.05 + 1.0 ) * probe.fwhm.c;
% 
%         tmp0 = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, [ probe.fwhm.r, probe.fwhm.c ] );
% 
%         tmp0 = tmp0 / max( abs( tmp0( : )));
% 
%         sol.probe.phi( :, :, pp ) = tmp0 .* exp( 0.0 * 2 * pi * 1i * tmp0 );
% 
%     end

else

    %=================================================
    % number of spatially coherent probe modes (SCPMs)
    %=================================================

    sol.probe.scpm.N = single( 5 ); 

    %========================================
    % Initialize to Gaussians of random sizes
    %========================================

    sol.probe.phi = zeros( [ sol.sz.sz, sol.probe.scpm.N ], 'single' );

    for pp = 1 : sol.probe.scpm.N

        probe.fwhm.r = single( 0.8e-6 / sol.csys.z2.dLy );
        probe.fwhm.c = single( 0.8e-6 / sol.csys.z2.dLx );

        probe.fwhm.r = (( 2 * rand - 1 ) * 0.05 + 1.0 ) * probe.fwhm.r;
        probe.fwhm.c = (( 2 * rand - 1 ) * 0.05 + 1.0 ) * probe.fwhm.c;

        tmp0 = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, [ probe.fwhm.r, probe.fwhm.c ] );

        tmp0 = tmp0 / max( abs( tmp0( : )));

        sol.probe.phi( :, :, pp ) = tmp0 .* exp( 0.0 * 2 * pi * 1i * tmp0 );

    end

    %========================
    % Guess of SCPM occupancy
    %========================

    sol.probe.scpm.occ = transpose( exp( +1 * 3.1 * ( 1 : sol.probe.scpm.N )));         % decaying exponential occupancy guess
    sol.probe.scpm.occ = sort( sol.probe.scpm.occ / norm( sol.probe.scpm.occ, 1 ));     % make sure the mode occupancy adds up to 1.0
    
    sol.probe.scpm.occ = single( sol.probe.scpm.occ );
    
    %==============
    % Orthogonalize
    %==============

    [ sol.probe.phi ] = orthog_modes_eigendecomp( sol.probe.phi );

    %========
    % Rescale
    %========

    [ sol.probe.phi ] = enforce_scpm_fro2TOT_photonocc( sol.probe.phi, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ );

end

sol.probe.phi    = single( sol.probe.phi );
sol.probe.scpm.N = single( sol.probe.scpm.N );

end