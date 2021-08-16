%
%{

clear; close all; testing_extract_replace_indices_noloops

%}

%==================================================================================================

clear; close all;
load /home/ash/Documents/Science/Matlab/Code/cdi/sim_ptycho2DTPA.mat;

[ sol, expt ] = set_defaults( sol, expt );  

sol.probe  = expt.probe;
sol.sample = expt.sample;
sol.spos   = expt.spos;
sol.sz     = expt.sz;

clearvars -except expt sol

sol.spos.updateorder = randperm( length( sol.spos.indxsubset ), round( 1.0 * length( sol.spos.indxsubset )));

%================================================
% what are the indices in the sample array for each scan position?

% [ ind ] = get_indices_2Dframes( expt );
% 
% ind   = gpuArray( ind );
% TFvec = gpuArray( expt.sample.TF( : ));
% 
% TFview3 = TFvec( ind );
% TFview4 = reshape( TFview3, [ expt.sz.sz, expt.spos.N ]);
% 
% TFview3g = gather( TFview3 );
% TFview4g = gather( TFview4 );

%================================================

sol.phi = zeros( [ sol.sz.sz, sol.probe.scpm.N, sol.spos.N ], 'single' );

for ss = 1 : sol.spos.N
    [ sol.phi( :, :, :, ss ), ~ ] = enforce_2DTPAsposview( expt.probe.P, expt.sample.TF, expt.sample.vs, expt.spos.rs( ss, : ), 'px' );
end
    

sol.spos.sampleindices = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs );

sol.probe.P = expt.probe.P;
sol.probe.scpm = expt.probe.scpm;

%==================================================================================================

tic

    % put back frames into 2D sample array
    
    [ TFm2 ] = DMupdate_2Dsample_loop( sol );
    
    
%     phi    = gpuArray( sol.phi );
%     probe  = gpuArray( sol.probe.P  );
%     updateorder = gpuArray( sol.spos.updateorder );
%     samsz  = gpuArray( sol.sample.sz.sz );
%     vs.r = gpuArray( sol.sample.vs.r );
%     vs.c = gpuArray( sol.sample.vs.c );
%     spos_rs = gpuArray( sol.spos.rs );
%     
%     [ TFm2 ] = DMupdate_2Dsample_loop( phi, probe, updateorder, samsz, vs, spos_rs );

toc
    


%================================================


    % put back frames into 2D sample array
    
    phi    = gpuArray( sol.phi );
    probe  = gpuArray( sol.probe.P  );
    ind    = gpuArray( sol.spos.sampleindices );
    sz     = gpuArray( sol.sz.sz );
    samsz  = gpuArray( sol.sample.sz.sz );
    spos_N = gpuArray( sol.spos.N );
    
tic

%     [ TFm ] = DMupdate_2Dsample_repmat( sol.phi, sol.probe.P, sol.spos.sampleindices, sol.sz.sz, sol.sample.sz.sz, sol.spos.N );
    [ TFm ] = DMupdate_2Dsample_repmat( phi, probe, ind, sz, samsz, spos_N );
    
toc

%==================================================================================================
    
    
    

%     % extract frames from 2D sample array
% %     [ ind ] = get_indices_2Dframes( sol );
%     [ ind ] = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs );
%     
% 
%     TFv0     = sol.sample.TF( : );
%     TFview_v = TFv0( ind );
%     TFview_m = reshape( TFview_v, [ expt.sz.sz, expt.spos.N ]);

    
    
    
    
    


    




    
    figure; imagesc( abs( sol.sample.TF ), [0, 1] )
    figure; imagesc( abs( TFm ), [0, 1]); title('repmat')
    figure; imagesc( abs( TFm2 ), [0, 1]); title('loop')
    
    
    norm( ( TFm - sol.sample.TF ) .* ( TFm ~=0 ), 'fro' )
    
    
    
    return
    
    

% t = 1 : ( 5 * 4 * 3 );
% t = reshape( t, 5, 4, 3 );



































































function [ sol, expt ] = set_defaults( sol, expt )

    warning('off','MATLAB:prnRenderer:opengl');
    
    rng( 'shuffle' )
 
    %================================================================
    %--------- Gaussian LPF for reconstruction resolution -----------
    %================================================================
        
    mu    = 0.5 * sol.sz.sz + 1;
    stdev = 0.80 * sol.sz.sz;
    tmp0  = make_2Dgaussian( sol.sz.sz, mu, stdev );
    
    sol.measLPF = 0 + 1 * fftshift( tmp0 );

    %============================================
    %---------- Fixed sample support ------------
    %============================================
    
%     sol.sample.support = 0 + 1 * make_rectangle( sol.sample.sz.sz, [ 0.43 * sol.sample.sz.r, 0.43 * sol.sample.sz.c ]);
%     % sol.probe.support = 0 + 1 * make_ellipsoid( sol.sample.sz.sz, [ 0.75 * sol.sample.sz.c, 0.75 * sol.sample.sz.r ]);
% 
%     [ sol.sample.support ] = 1 + 0 * lpf_gauss( sol.sample.support, 222.01 * sol.sample.sz.sz );
 
    %{
    
    figure; imagesc( sol.sample.support .* angle( sol.sample.TF ))
    figure; imagesc( sol.sample.support )
    
    %}
    

    %======================================================================
    %---------- Sample magnitude/phase scaling + modifications ------------
    %======================================================================
    
    %{
    
    abs_TF = abs(  sol.sample.TF );
    abs_TF = abs_TF / max( abs_TF( : ));
    
    phs_TF = angle(  sol.sample.TF );
    phs_TF = phs_TF - min( phs_TF( : ));
    phs_TF = phs_TF / max( phs_TF( : ));
    
    as = 0.2;
    ps = 0.5;
    
    tmp0 = ( as * abs_TF + ( 1 - as ) * phs_TF ) .* exp( 1i * 2 * pi * (  ps * abs_TF + ( 1 - ps ) * phs_TF ));
    
    sol.sample.TF = tmp0;

    %}



    %{

    abs_TF = abs( sol.sample.TF );
    phs_TF = angle( sol.sample.TF );

    tmp0 = abs_TF .* exp( 1i * phs_TF + 1i * 0.0 );

    figure; 
    
    ax1 = subplot(131); imagesc( abs( tmp0 ));   daspect([1 1 1]); colormap( ax1, expt.cm.blkgrn )
    ax2 = subplot(132); imagesc( angle( tmp0 )); daspect([1 1 1]); colormap( ax2, expt.cm.blj )
    subplot(133);       imagescHSV( tmp0 );      daspect([1 1 1]);

    sol.sample.TF = tmp0;

    %}

    %============================================
    %---------- Fixed probe support -------------
    %============================================

    sol.probe.support = 0 + 1 * make_rectangle( sol.sz.sz, [ 0.9 * sol.sz.r, 0.9 * sol.sz.c ]);
    % sol.probe.support = 0 + 1 * make_ellipsoid( sol.sz.sz, [ 0.75 * sol.sz.c, 0.75 * sol.sz.r ]);

    [ sol.probe.support ] = 1 + 0 * lpf_gauss( sol.probe.support, 0.03 * sol.sz.sz );

    %======================================================
    %------ When to start sample/probe/spos updates -------
    %======================================================
    
    sol.it.spos_start  = 1e99;
    sol.it.spos_update = 1e99;
    
    sol.spos.suofai  = logical( 0 );                 % same update order for all iterations
    sol.spos.rand_ss = 1/2;                         % 

    sol.it.sample_start        = 0;
    sol.it.sample_update       = 1;
    sol.it.sample_sparsity     = 1;
    sol.it.sample_mag_phs_ineq = 1;
    
    sol.it.probe_start   = 0;
    sol.it.probe_update  = 1;
    sol.it.probe_orthog  = 1;
    sol.it.probe_scaling = 1;
    sol.it.probe_maxvals = 1;
    sol.it.probe_support = 1;

    sol.it.print_img_results = 10;
    sol.it.collect_metrics   = 50;

    %============================================
    %-------- PR projection algo params ---------
    %============================================
     
    sol.RAAR.beta = 0.5;

    %============================================
    %----- Sparsity constraints for sample ------
    %============================================
    
    sol.sparse.s2DFDxy = setup_edgedetect_forwarddiff2Dxy( [], sol.sample.sz.sz );
    
    sol.sparse.type = 'sparse2DFDxy';
    
    sol.sparse.type_sparse2DFDxy = { 'aniso_abs', ...                  % 1
                                     'aniso_phs', ...                  % 2
                                     'iso_abs', ...                    % 3
                                     'iso_phs', ...                    % 4
                                     'iso_abs_iso_phs', ...            % 5
                                     'aniso_abs_aniso_phs', ...        % 6
                                     'aniso_re_aniso_im', ...          % 7
                                     'iso_re_iso_im' };                % 8

    %============================================
    %-------------- Sample scaling --------------
    %============================================
    
    sol.sample.phsL = 0;
    sol.sample.phsH = +1.9 * pi;

    sol.sample.absL = 0.0;
    sol.sample.absH = 1.0;
    
    %============================================
    %-------------- Probe scaling ---------------
    %============================================
    
    %%{

    sol.probe.P = orthog_modes_eigendecomp( sol.probe.P );
    
    %========================

    sol.probe.scpm.mean_meas_fro2TOT = mean( expt.meas.SI_sumD2 );
%     sol.probe.scpm.fro2TOT = 4.20 * sol.probe.scpm.mean_meas_fro2TOT;
%     sol.probe.scpm.fro2TOT = 7e4 ^ 2;
    sol.probe.scpm.fro2TOT = 1400 ^ 2; 

    %========================
    
    sol.probe.scpm.max = 10 * [ 10, 20, 50 ];

    %========================
    
    sol.probe.scpm.occ = [ ];
%     sol.probe.scpm.occ = [ 0.05, 0.15, 0.80 ];
%     sol.probe.scpm.occ = [ 0.05, 0.10, 0.85 ];
%     sol.probe.scpm.occ = [ 0.03, 0.07, 0.90 ];
    sol.probe.scpm.occ = sort( sol.probe.scpm.occ / norm( sol.probe.scpm.occ, 1 ));     % make sure the mode occupancy adds up to 1.0
    
    sol.probe.P = enforce_scpm_fro2TOT_photonocc( sol.probe.P, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ );
    
    %}
    
    %============================================
    %---------- Scan positions init -------------
    %============================================
    
    sol.spos.update.dpx        = 1/4;  % for computing finite differences wrt scan position
    sol.spos.update.maxpx      = 4 * sol.spos.update.dpx;
    sol.spos.update.linesearch = transpose( 0.0 : sol.spos.update.dpx : sol.spos.update.maxpx );
    
%     sol.spos.update.shifttype = 'px';
    sol.spos.update.shifttype = 'subpx';

%     sol.spos.update.grad = 'all';
    sol.spos.update.grad = 'indiv';

    % for this data set, we use 0.5 um steps...only allow us to go so far from initial spos we start at
    sol.spos.update.maxcorrectr = 1.00e-6 / expt.csys.z2.dLv;
    sol.spos.update.maxcorrectc = 1.00e-6 / expt.csys.z2.dLh;
    
    sol.spos.indxsubset = expt.spos.indxsubset;                                                 % the scan positions we update the exitwaves/sample/probes over
    sol.spos.shifttype = 'px';                                                                  % When extracting part of the sample in the current scan position, use this type of pixel shifting  

    %================================================================
    %------ Center the probe function using center of mass ----------
    %================================================================
    
    %{

    figure; imagesc(abs(sol.probe.P( :, :, end )))

    [ com ] = centerOfMass( abs( sol.probe.P( :, :, end )));

    for pp = 1 : sol.probe.scpm.N

        sol.probe.P( :, :, pp ) = circ shift( sol.probe.P( :, :, pp ), -1 * round( com - 0.5 * sol.sz.sz - 1 ));

    end

    figure; imagesc(abs(sol.probe.P( :, :, end )))


    [~, Ic ] = max( sum( abs(sol.probe.P( :, :, end )), 1 ));
    [~, Ir ] = max( sum( abs(sol.probe.P( :, :, end )), 2 ));

    sol.probe.P( :, :, pp ) = nocircshift2D( sol.probe.P( :, :, pp ), -1 * round( [ Ir, Ic ] - 0.5 * sol.sz.sz - 1 ));

    figure; imagesc(abs(sol.probe.P( :, :, end )))

    %}


    %============================================
    %--------- Compute initial exit waves  ------
    %============================================

    
    
%     [ sol.phi ]       = ERupdate_2DTPAexwv_spos_meas( sol, expt );  
%     [ sol.sample.TF ] = ePIEupdate_sample( sol.probe.P, sol.sample.TF, sol.phi, sol );
%     [ sol.probe.P ]   = DMupdate_probemodes( sol.sample.TF, sol.phi, sol ); 
%     
%     sol.sample.TF = modulus_limits_project( sol.sample.TF, [ sol.sample.absL, sol.sample.absH ] );
%     % sol.sample.TF = modulus_limits_scale( sol.sample.TF, [ sol.sample.absL, sol.sample.absH ] );
%     sol.sample.TF = phase_limits_project( sol.sample.TF, [ sol.sample.phsL, sol.sample.phsH ] );   
            


    sol.phi = zeros( [ sol.sz.sz, sol.probe.scpm.N, sol.spos.N ], 'single' );

    for ss = 1 : sol.spos.N

        rs = +sol.spos.rs( ss, : );
    %     rs = -sol.spos.rs( ss, : );

        [ sol.phi( :, :, :, ss ), ~ ] = enforce_2DTPAsposview( sol.probe.P, sol.sample.TF, sol.sample.vs, +1 * rs, sol.spos.shifttype );

    end

    
    
    

end


    
    