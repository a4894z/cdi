function [ sol, expt ] = runsolver_ptycho2DTPA

    %================================================================================================================================================

    %{

    cd /net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi

    %========

    codelocation = '~/Documents/Science/Matlab/Code/cdi/'; cd( codelocation ); restoredefaultpath; addpath( genpath( pwd ));

    %========

    clear; close all; 

    restoredefaultpath

    % codelocation =  '~/Documents/Science/Matlab/Code/cdi/';

    codelocation =  '/net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi';
    cd( codelocation );

    addpath( genpath( pwd ));   

    clearvars -except expt sol

    %========

    % clear; close all; 
    
    clear; close all; [ sol, expt] = runsolver_ptycho2DTPA
    
    %========

    %}

    %================================================================================================================================================
    
    rng( 'shuffle' )
    
    restoredefaultpath; 
    addpath( genpath( pwd ));
    clearvars -except expt sol
    
    %================================================================================================================================================
    
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/noise/sim_ptycho2DTPA.mat';    % WITH    NOISE
    data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA.mat';   % WITHOUT NOISE

%     data_path = [ pwd, '/sim_ptycho2DTPA.mat' ];

    %=========
    
    load( data_path, 'sol', 'expt' );  

    expt.paths.rsdata = [ pwd, '/sim_ptycho2DTPA.mat' ];

    clearvars -except expt sol

    %================================================================================================================================================

    % Initialize default parameters for use in phase retrieval
    [ sol, expt ] = runsolver_ptycho2DTPA_config( sol, expt );    
    
    %=========
    
    sol.rPIE_alpha = single( 0.10 );
    
    %=========
    
%     sol.sample.T = rand( sol.sample.sz.sz, 'single' ) + 1i * rand( sol.sample.sz.sz, 'single' );
%     sol.sample.T = sol.sample.T / max( abs( sol.sample.T( : )));
    
    %=========
    % GPU Init
    %=========
    
    sol.use_gpu = true; 

    sol.gpu_id = 1; 
 
    if sol.use_gpu == true, reset( gpuDevice( sol.gpu_id )); end

    %================================
    % Stochastic minibatch parameters
    %================================

    sol.spos.rand_spos_subset_pct = 0.10;      
    
    sol.spos.rand_spos_subset_pct = single( sol.spos.rand_spos_subset_pct );  
    if sol.spos.rand_spos_subset_pct > 1.00, sol.spos.rand_spos_subset_pct = single( 1.00 ); end

    %======================================================
    % epoch, iteration, and iteration repeat specifications
    %======================================================
    
    N_pauseandsave              = single( 1 );
    
    N_epochs                    = single( 5000 );

    sol.it.probe_orthog         = single( 25 );
    
    sol.it.metrics_and_plotting = single( 100 );
    sol.print_img               = logical( 0 );             %#ok<LOGL>

    %============
    % CHEAT CODES
    %============

    % sol.sample.T      = expt.sample.T;
    % sol.probe.phi     = expt.probe.phi;
    % sol.probe.scpm     = expt.phi.scpm;
    % sol.probe.scpm.max = 1e3; 

    %================================================================================================================================================

    for ii = 1 : N_pauseandsave

        %==========
        % minibatch
        %==========

        [ sol, expt ] = ptycho2DTPA_runGPU_stochminibatchgrad( sol, expt, N_epochs );    

        %==========
        % fullbatch
        %==========

%         [ sol, expt ] = ptycho2DTPA_runGPU_stochcoordgrad( sol, expt, N_epochs );         
%         [ sol, expt ] = ptycho2DTPA_runGPU_totalgrad( sol, expt, N_epochs );     

        %========

        % Clean up   
        [ sol, expt ] = ptycho2DTPA_cleanup( sol, expt );     

        % Save results 
        expt = ptycho2DTPA_saveresults( sol, expt, ii );    

        % clear GPU memory
        if sol.use_gpu == true, reset( gpuDevice( sol.gpu_id )); end

    end

end

%====================================================================================================================================================

function [ sol, expt ] = runsolver_ptycho2DTPA_config( sol, expt )     

    warning('off','MATLAB:prnRenderer:opengl');
    
    rng( 'shuffle' )
             
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % projection algorithm parameter(s)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sol.RAAR_beta  = single( 0.5 );
    sol.rPIE_alpha = single( 0.05 );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gaussian LPF for reconstruction resolution 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    mu    = single( 0.5 * sol.sz.sz + 1 );
    stdev = single( 0.90 * sol.sz.sz );
    tmp0  = make_2Dgaussian( sol.sz.sz, mu, stdev );
    
    sol.measLPF = single( 1 + 0 * fftshift( tmp0 ));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fliplr/flipud/rot90() of measurements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     expt.meas.D = nocircshift3D( expt.meas.D, [ 0, -1, 0 ] );

%     expt.meas.D    = fliplr( expt.meas.D );%     expt.meas.D    = rot90( expt.meas.D, -1 );
%     expt.meas.D    = fliplr( expt.meas.D );
% 
%     expt.meas.Deq0 = ( expt.meas.D == 0 );

%     expt.meas.D    = flipud( expt.meas.D );



%     % FOR L0013   
% %     expt.meas.D = rot90( expt.meas.D, -1 ); 
%     expt.meas.D = flipud( expt.meas.D );  
    
%     % FOR L0011   
%     expt.meas.D = rot90( expt.meas.D, -1 ); 
%     expt.meas.D = fliplr( expt.meas.D );          
% 
%     expt.meas.Deq0 = ( expt.meas.D == 0 );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % When to start sample/probe/spos updates 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sol.it.spos_start  = single( 1e99 );
    sol.it.spos_update = single( 1e99 );

    %========
    
    sol.it.metrics_and_plotting = single( 5 );

    sol.it.sample_start        = single( 0 );
    sol.it.sample_update       = single( 1 );
    sol.it.sample_mag_phs_ineq = single( 1 ); 
    sol.it.sample_sparsity     = single( 1e99 );
    
    sol.it.probe_start     = single( 0 );      
    sol.it.probe_update    = single( 1 );
    sol.it.probe_scaling   = single( 1 );
    sol.it.probe_support   = single( 1 );
    sol.it.probe_centering = single( 1e99 );
    sol.it.probe_orthog    = single( 50 );
    sol.it.probe_maxvals   = single( 1e99 );

    %%%%%%%%%%%%%%%%%%%%%%%
    % PROBE INITIALIZATIONS 
    %%%%%%%%%%%%%%%%%%%%%%%
    
    %========================================
    % Tweak spatially incoherent probe modes?
    %========================================

%     Pmodes_old  = sol.probe.phi;
%     sol.probe.phi = [];
%     
%     sol.probe.phi( :, :, 5 ) = Pmodes_old( :, :, 3 );
%     sol.probe.phi( :, :, 4 ) = Pmodes_old( :, :, 2 );
%     sol.probe.phi( :, :, 3 ) = Pmodes_old( :, :, 1 );
%   
%     sol.probe.phi( :, :, 2 ) = single( 0.02 + 0 * sol.probe.phi( :, :, 2 ));
%     sol.probe.phi( :, :, 1 ) = single( 0.01 + 0 * sol.probe.phi( :, :, 1 ));
%     
%     sol.probe.scpm.N = size( sol.probe.phi, 3 );
%     
% %     sol.probe.scpm.occ = [ 0.02, 0.03, 0.95 ];
% %     sol.probe.scpm.occ = sort( sol.probe.scpm.occ / norm( sol.probe.scpm.occ, 1 ));     % make sure the mode occupancy adds up to 1.0
%      
%  %     sol.probe.phi = enforce_scpm_fro2TOT_photonocc( sol.probe.phi, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ );
%      sol.probe.phi = orthog_modes_eigendecomp( sol.probe.phi );
     
    %=================================================
    % Probe initializations: Max abs probe mode values
    %=================================================
    
%     sol.probe.scpm.max = [ 400, 400, 400 ];
    sol.probe.scpm.max = single( [ ] );

    try sol.probe.scpm.max = reshape( sol.probe.scpm.max, [ 1, 1, sol.probe.scpm.N ] ); catch, end
    
    %============================================
    % Probe initializations: Probe mode occupancy
    %============================================
    
    sol.probe.scpm.occ = single( [ ] );
%     sol.probe.scpm.occ = single( [ 0.10, 0.20, 0.70 ] ); 
%     sol.probe.scpm.occ = single( [ 0.10, 0.25, 0.65 ] ); 
%     sol.probe.scpm.occ = single( [ 0.05, 0.15, 0.80 ] );
%     sol.probe.scpm.occ = single( [ 0.05, 0.10, 0.85 ] );
%     sol.probe.scpm.occ = single( [ 0.03, 0.07, 0.90 ] );
    
    sol.probe.scpm.occ = sort( sol.probe.scpm.occ / norm( sol.probe.scpm.occ, 1 ));     % make sure the mode occupancy adds up to 1.0
      
%     sol.probe.phi = orthog_modes_eigendecomp( sol.probe.phi );

    %=========================
    % probe mode total scaling
    %=========================
    
%     sol.probe.scpm.fro2TOT = single( [];
    sol.probe.scpm.fro2TOT = single( 1.75 * mean( expt.meas.SI_sumD2 ));
%     sol.probe.scpm.fro2TOT = single( 1.05 * ( 3e3 ^ 2 ));   
%     sol.probe.scpm.fro2TOT = single( 1.10 * expt.phi.scpm.fro2TOT);   
    
    sol.probe.phi = enforce_scpm_fro2TOT_photonocc( sol.probe.phi, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ );
    
    %====================
    % Fixed probe support 
    %====================
    
    sol.probe.support = make_rectangle( sol.sz.sz, [ 0.8 * sol.sz.r, 0.8 * sol.sz.c ]);
    % sol.probe.support = make_ellipsoid( sol.sz.sz, [ 0.75 * sol.sz.c, 0.75 * sol.sz.r ]);
    
%     sol.probe.support = lpf_gauss( sol.probe.support, 90.03 * sol.sz.sz );
    
    sol.probe.support( abs( sol.probe.support ) < 1e-3 ) = 0;
    
    sol.probe.support = 0 + 1 * sol.probe.support;
    
%     sol.probe.phi = sol.probe.phi .* sol.probe.support;
    
    %=========================
    
    sol.probe.scpm.N       = single( sol.probe.scpm.N );
    sol.probe.scpm.fro2TOT = single( sol.probe.scpm.fro2TOT );
    sol.probe.scpm.max     = single( sol.probe.scpm.max );
    sol.probe.scpm.occ     = single( sol.probe.scpm.occ );
    
    %=========================
    % Shrinkwrap probe support 
    %=========================
    
    sol.swparams.blurx     = single( 0.02 ); 
    sol.swparams.blury     = single( 0.02 );
    sol.swparams.sparselvl = single( 0.60 ); 
                        
    %===========================
    % Center the probe intensity
    %===========================
    
    %{

    figure; imagesc(abs(sol.probe.phi( :, :, end )))

    [ com ] = centerOfMass( abs( sol.probe.phi( :, :, end )));

    for pp = 1 : sol.probe.scpm.N

        sol.probe.phi( :, :, pp ) = circ shift( sol.probe.phi( :, :, pp ), -1 * round( com - 0.5 * sol.sz.sz - 1 ));

    end

    figure; imagesc(abs(sol.probe.phi( :, :, end )))


    [~, Ic ] = max( sum( abs(sol.probe.phi( :, :, end )), 1 ));
    [~, Ir ] = max( sum( abs(sol.probe.phi( :, :, end )), 2 ));

    sol.probe.phi( :, :, pp ) = nocircshift2D( sol.probe.phi( :, :, pp ), -1 * round( [ Ir, Ic ] - 0.5 * sol.sz.sz - 1 ));

    figure; imagesc(abs(sol.probe.phi( :, :, end )))

    %}

    %%%%%%%%%%%%%%%%%%%%%%%
    % Sample initialization 
    %%%%%%%%%%%%%%%%%%%%%%%

    %================================
    % Fixed sample support/shrinkwrap
 
%     sol.sample.support = 0 + 1 * make_rectangle( sol.sample.sz.sz, [ 0.43 * sol.sample.sz.r, 0.43 * sol.sample.sz.c ]);
%     % sol.probe.support = 0 + 1 * make_ellipsoid( sol.sample.sz.sz, [ 0.75 * sol.sample.sz.c, 0.75 * sol.sample.sz.r ]);
% 
%     [ sol.sample.support ] = 1 + 0 * lpf_gauss( sol.sample.support, 222.01 * sol.sample.sz.sz );
 
    %{
    
    figure; imagesc( sol.sample.support .* angle( sol.sample.T ))
    figure; imagesc( sol.sample.support )
    
    %}

    %===============
    % Sample scaling 
 
    sol.sample.phsL = single( -0.50 * pi );
    sol.sample.phsH = single( +0.50 * pi );

    sol.sample.absL = single( 0.0 );
    sol.sample.absH = single( 1.1 );
   
    %================================
    % Sparsity constraints for sample 

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
   
    % only use sample sparsity on subregion of sample image
    sol.sparse.support = ones( sol.sample.sz.sz, 'logical' ); 

    sol.sparse.pct         = single( 0.20 );
    sol.sparse.lvl         = round( sol.sample.sz.rc * sol.sparse.pct );
    sol.sparse.threshtype  = 's';
    sol.sparse.threshname  = sol.sparse.type_sparse2DFDxy{ 1 };

    %====================================================
    % Sample magnitude/phase scaling + misc modifications
 
    %{
    
    abs_TF = abs(  sol.sample.T );
    abs_TF = abs_TF / max( abs_TF( : ));
    
    phs_TF = angle(  sol.sample.T );
    phs_TF = phs_TF - min( phs_TF( : ));
    phs_TF = phs_TF / max( phs_TF( : ));
    
    as = 0.2;
    ps = 0.5;
    
    tmp0 = ( as * abs_TF + ( 1 - as ) * phs_TF ) .* exp( 1i * 2 * pi * (  ps * abs_TF + ( 1 - ps ) * phs_TF ));
    
    sol.sample.T = tmp0;

    %}
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Scan Position Correction
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sol.spos.update.dpx        = single( 1/4 );  % for computing finite differences wrt scan position
    sol.spos.update.maxpx      = single( 4 * sol.spos.update.dpx );
    sol.spos.update.linesearch = single( transpose( 0.0 : sol.spos.update.dpx : sol.spos.update.maxpx ));
    
%     sol.spos.update.shifttype = 'px';
    sol.spos.update.shifttype = 'subpx';

%     sol.spos.update.grad = 'all';
    sol.spos.update.grad = 'indiv';

    % for this data set, we use 0.5 um steps...only allow us to go so far from initial spos we start at
    sol.spos.update.maxcorrectr = single( 1.00e-6 / expt.csys.z2.dLy );
    sol.spos.update.maxcorrectc = single( 1.00e-6 / expt.csys.z2.dLx );
    
    sol.spos.indx = single( expt.spos.indx );       % the scan positions indices we update the exitwaves/sample/probes over
    sol.spos.shifttype = 'px';                      % When extracting part of the sample in the current scan position, use this type of pixel shifting  

end

