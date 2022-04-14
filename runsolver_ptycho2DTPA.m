function [ sol, expt ] = runsolver_ptycho2DTPA

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %{

    cd /net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi

    %========

    clear; close all; 

    restoredefaultpath

    % codelocation =  '~/Documents/Science/Matlab/Code/cdi/';
    codelocation =  '/net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi';
    
    cd( codelocation );

    addpath( genpath( pwd ));   

    clearvars -except expt sol

    %========

    clear; close all; [ sol, expt] = runsolver_ptycho2DTPA;
    
    %========
    
    /net/s8iddata/export/8-id-ECA/Analysis/reduced_data/2021-3/zjiang20211214
    /net/s8iddata/export/8-id-i/2020-2/zjiang202007/reduced_data

    %}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %#ok<*LOGL>
    warning('off','MATLAB:prnRenderer:opengl');
    
    rng( 'shuffle' )
    
    %==============================================================
    % Set paths to code location and add relevant folders and files
    %==============================================================
    
    restoredefaultpath; 
    addpath( genpath( pwd ));
%     addpath( genpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/' ));
%     addpath( genpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/' ));
    clearvars -except expt sol
    
    %============================
    % Set paths to load data from
    %============================
    
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/noise/sim_ptycho2DTPA.mat';      % WITH NOISE, WITHOUT BG RM
    data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/noise_rmbg/sim_ptycho2DTPA.mat'; % WITH NOISE, WITH BG RM
    
%     data_path = [ pwd, '/sim_ptycho2DTPA.mat' ];

    %=========
    
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p1degree_1088x2560.mat';  
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p2degree_1088x2560.mat'; 
    
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p3degree_1088x2560.mat';  
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p4degree_1088x2560.mat'; 
    
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p5degree_1088x2560.mat';  
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p6degree_1088x2560.mat';  
    
    %========================
    % Load the processed data 
    %========================
    
    load( data_path, 'sol', 'expt' );  

%     expt.paths.rsdata = data_path;
    expt.paths.rsdata = [ './', expt.paths.data_mat_name, '.mat' ];
%     expt.paths.rsdata = './pecoCSSIsimulatedmultislicedata.mat';
    
    clearvars -except expt sol
    
    %===============================================================
    % Initialize/reset default parameters for use in phase retrieval
    %===============================================================
    
    sol.init.reinit_sampleTF = logical( 0 );
     
    sol.init.reinit_scpm_max = logical( 0 );   % reinit the scpm max magnitudes
    sol.init.reinit_scpm_occ = logical( 0 );   % reinit the scpm occupancies
    sol.init.reinit_scpm_phi = logical( 0 );   % reinit the scpm shapes

    [ sol, expt ] = runsolver_ptycho2DTPA_initialize_phaseretrieval_parameters( sol, expt );   
    
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

    %======================================
    % !!!! CHEAT CODES FOR SIMULATIONS !!!!
    %======================================

%     sol.sample.T = expt.sample.T;
%     
%     sol.probe.phi      = expt.probe.phi;
%     sol.probe.scpm.occ = expt.probe.scpm.occ;
%     sol.probe.scpm.N   = expt.probe.scpm.N;
%     sol.probe.scpm.max = max( reshape( abs( expt.probe.phi ), [ expt.sz.rc, expt.probe.scpm.N ] ), [], 1 );
    
%     sol.spos.rs    = expt.spos.rs;
%     sol.spos.indx  = expt.spos.indx;
    
    %================================================================================================================================================
    
    N_pauseandsave = single( 1 );
    N_epochs       = single( 500 );
    
    for ii = 1 : N_pauseandsave

        %=============
        % minibatch GD
        %=============

        [ sol, expt ] = ptycho2DTPA_runGPU_stochminibatchgrad( sol, expt, N_epochs );    

        %==============================
        % fullbatch block stochastic GD
        %==============================

%         [ sol, expt ] = ptycho2DTPA_runGPU_stochcoordgrad( sol, expt, N_epochs );    ( !!!!!!!!!! DON'T USE, NEED TO UPDATE !!!!!!!!!! )

        %===================
        % fullbatch total GD 
        %===================

%         [ sol, expt ] = ptycho2DTPA_runGPU_totalgrad( sol, expt, N_epochs );       ( !!!!!!!!!! DON'T USE, NEED TO UPDATE !!!!!!!!!! )

        %==========================
        % bookkeeping and logistics
        %==========================
        
        % Clean up   
        [ sol, expt ] = ptycho2DTPA_cleanup( sol, expt );     

        % Save results 
        expt = ptycho2DTPA_saveresults( sol, expt, ii );    

        % Clear/reset GPU memory
        if sol.use_gpu == true, reset( gpuDevice( sol.gpu_id )); end

    end

end

%====================================================================================================================================================

function [ sol, expt ] = runsolver_ptycho2DTPA_initialize_phaseretrieval_parameters( sol, expt )     
    
    %=====================================
    % rPIE adaptive step length parameters
    %=====================================
    
    sol.rPIE_alpha_phi = single( 1e-0 );   % relatively insensitive to this?
    sol.rPIE_alpha_T   = single( 1e-0 );   
    
    %========================================
    % When to start sample/probe/spos updates 
    %========================================
    
    sol.it.spos_start  = single( 1e99 );
    sol.it.spos_update = single( 1e99 );

    sol.it.sample_start    = single( 0 );
    sol.it.sample_update   = single( 1 );
    sol.it.sample_mag_ineq = single( 1 ); 
    sol.it.sample_phs_ineq = single( 1e99 ); 
    sol.it.sample_sparsity = single( 5 );

    sol.it.probe_start        = single( 5 );      
    sol.it.probe_update       = single( 1 );
    sol.it.probe_scaling      = single( 1 );
    sol.it.probe_support      = single( 1 );
    sol.it.probe_centering    = single( 1 );
    sol.it.probe_orthog       = single( 50 );   % 25 for rPIE paper
    sol.it.probe_maxvals      = single( 1 );
    sol.it.probe_smoothing    = single( 5e99 );
    sol.it.probe_phase_sparse = single( 1e99 );
    
    sol.it.collect_metrics   = single( 10 );     % when we want to collect performance metrics
    sol.it.mkimg_meas_metric = single( 50 ); 
    sol.it.mkimg_sample_SCPM = single( 50 ); 
    
    %=============================================
    % Performance metric and epoch update counters
    %=============================================
    
    if ~isfield( sol.it, 'mtot' ),  sol.it.mtot  = 1; end
    if ~isfield( sol.it, 'metr' ),  sol.it.metr  = 1; end
    if ~isfield( sol.it, 'epoch' ), sol.it.epoch = 1; end

%     sol.it.metr  = 1;
%     sol.it.epoch = 1;
%     sol.it.exwv  = 1;

    % for metrics plotting
    sol.metrics.metrics_ylim      = [ 10^-3, 10^0 ];
    sol.metrics.grad_metrics_ylim = [ 10^-4, 10^0 ];
    
%     sol.metrics.legend_loc = 'northeast';
    sol.metrics.legend_loc = 'southwest';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCAN POSITION INITIALIZATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %======================================
    % Scan Position Shearing Transformation
    %======================================
    
%     rs = sol.spos.rs;
    
%     s_x = 0.9;
%     sol.spos.shear_x = [ 1, 0; s_x, 1 ];
%     
%     sol.spos.rs = transpose( sol.spos.shear_x * transpose( sol.spos.rs ));
%     
%     %========

%     sol.sample.T = padarray( sol.sample.T, [ 0, 1 ], 1, 'post' );
%     
%     sol.sample.sz.sz      = size( sol.sample.T );
%     sol.sample.sz.r       = sol.sample.sz.sz( 1 );
%     sol.sample.sz.c       = sol.sample.sz.sz( 2 );
%     sol.sample.sz.rc      = sol.sample.sz.r * sol.sample.sz.c;
%     sol.sample.sz.sqrt_rc = sqrt( sol.sample.sz.rc );
    
%     % determine proper view region FOV ( based on probe size ) for when computing the exit wave views at a particular scan position
%     sol.sample.vs.r = single( round( ( 0.5 * ( sol.sample.sz.r - sol.sz.r ) + 1 ) : ( 0.5 * ( sol.sample.sz.r + sol.sz.r ))));
%     sol.sample.vs.c = single( round( ( 0.5 * ( sol.sample.sz.c - sol.sz.c ) + 1 ) : ( 0.5 * ( sol.sample.sz.c + sol.sz.c ))));

%     figure; 
%     plot_2Dscan_positions( sol.spos.rs, [], sol.spos.rs, [] )
%     set( gca, 'xdir', 'reverse' )
%     set( gca, 'ydir', 'normal' )
% %     daspect([1 .1 1])
%     xlabel('xh, lab frame'); 
%     ylabel('yv, lab frame');
%     %title('positions for scanning the probe on the sample')

    %================================================
    % Gradient Descent Based Scan Position Correction
    %================================================
    
%     sol.spos.update.dpx        = single( 1/4 );  % for computing finite differences wrt scan position
%     sol.spos.update.maxpx      = single( 4 * sol.spos.update.dpx );
%     sol.spos.update.linesearch = single( transpose( 0.0 : sol.spos.update.dpx : sol.spos.update.maxpx ));
%     
% %     sol.spos.update.shifttype = 'px';
%     sol.spos.update.shifttype = 'subpx';
% 
% %     sol.spos.update.grad = 'all';
%     sol.spos.update.grad = 'indiv';
% 
%     % for this data set, we use 0.5 um steps...only allow us to go so far from initial spos we start at
%     sol.spos.update.maxcorrectr = single( 1.00e-6 / expt.csys.z2.dLy );
%     sol.spos.update.maxcorrectc = single( 1.00e-6 / expt.csys.z2.dLx );
    
    %=================================================================================================
    % When extracting part of the sample in the current scan position, use this type of pixel shifting  
    %=================================================================================================
    
    sol.spos.shifttype = 'px';   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXITWAVE/MEASUREMENT INITIALIZATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %=====================================================
    % Poisson or Gaussian Noise Model for Exitwave Updates
    %=====================================================
   
%     sol.exwv_noisemodel = 'gaussian';
    
    sol.exwv_noisemodel = 'poisson';
    sol.it.exwv_poisson_exactsearch = 1;

    %========================================================================================================
    % when monitoring the Poisson cost function metric, we need a scaling offset so that it's always positive
    %========================================================================================================
    
    if ~isfield( sol, 'metrics' ), sol.metrics = struct; end
        
    if ~isfield( sol.metrics, 'poisson_offset' )
    
        I_m                        = expt.meas.D .^ 2;
        log_I_m                    = log( I_m );
        log_I_m( isinf( log_I_m )) = 0;

        sol.metrics.poisson_offset = ( I_m - I_m .* log_I_m );
        sol.metrics.poisson_offset = sum( sol.metrics.poisson_offset(:) ) / size( I_m, 3 );

        clear( 'I_m', 'log_I_m' )

    end
    
    %===============================================================
    % Gaussian LPF for suppressing noise at high spatial frequencies
    %===============================================================
        
%     sol.measLPF = make_2Dgaussian( sol.sz.sz, single( 0.50 * sol.sz.sz + 1 ), single( 0.5 * sol.sz.sz ) );
    
    sol.measLPF = ones( sol.sz.sz, 'single' );
    
    %==========================================
    % fliplr/flipud/rot90()/etc of measurements
    %==========================================

%     expt.meas.D = flip( expt.meas.D, 2 );
%     expt.meas.D = flip( expt.meas.D, 1 );
%     expt.meas.D = nocircshift3D( expt.meas.D, [ 0, +2, 0 ] );

%     tmp0 = sol.spos.rs;
%     sol.spos.rs( :, 1 ) = tmp0( :, 2 );
%     sol.spos.rs( :, 2 ) = tmp0( :, 1 );
%     
%     clear( 'tmp0' )

%     expt.meas.D = nocircshift3D( fftshift( fftshift( expt.meas.D, 1 ), 2 ), [ +4, 0, 0 ] );
%     expt.meas.D = fftshift( fftshift( expt.meas.D, 1 ), 2 );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBE INITIALIZATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %============================
    % Gaussian smoothing of SCPMs
    %============================
    
    sol.probe.gaussian_blur_std_x = 0.05 * sol.sz.c;
    sol.probe.gaussian_blur_std_y = 0.05 * sol.sz.r;
    
    sol.probe.gaussian_lpf = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, [ sol.probe.gaussian_blur_std_y, sol.probe.gaussian_blur_std_x ] );
    
    sol.probe.gaussian_lpf = fftshift( sol.probe.gaussian_lpf );
     
    %================================
    % Max abs probe mode pixel values
    %================================
    
    if ~isfield( sol.probe.scpm, 'max' ), sol.probe.scpm.max = single( [ ] ); end
    
    %========
    
    if isempty( sol.probe.scpm.max ) || sol.init.reinit_scpm_max
    
        sol.probe.scpm.max = single( exp( +0.1 * ( 1 : sol.probe.scpm.N )));
        
        sol.probe.scpm.max = single( 200 * sol.probe.scpm.max / max( sol.probe.scpm.max( : )));
    
    end
    
    % try to reshape, do nothing if the max val array is empty
    try sol.probe.scpm.max = single( reshape( sol.probe.scpm.max, [ 1, 1, sol.probe.scpm.N ] )); catch, end
    
    %=====================
    % Probe mode occupancy
    %=====================
    
    if ~isfield( sol.probe.scpm, 'occ' ), sol.probe.scpm.occ = single( [ ] ); end
    
    %========
    
    if isempty( sol.probe.scpm.occ ) || sol.init.reinit_scpm_occ 
        
%         sol.probe.scpm.occ = transpose( [ 0.25, 0.75 ]);              
        sol.probe.scpm.occ = transpose( exp( +1 * 0.6 * ( 1 : sol.probe.scpm.N )));                 % decaying exponential occupancy guess
        
        sol.probe.scpm.occ = single( sort( sol.probe.scpm.occ / norm( sol.probe.scpm.occ, 1 )));     % make sure the mode occupancy adds up to 1.0
        
    end

    %========================
    % probe intensity scaling
    %========================
    
    sol.probe.scpm.fro2TOT = single( [] );
    
    sol.probe.scpm.fro2TOT = single( 1.75 * mean( squeeze( sum( sum( expt.meas.D .^ 2, 1 ), 2 ))));   % 1.75 for rPIE paper
%     sol.probe.scpm.fro2TOT = single( 0.9 * expt.probe.scpm.fro2TOT );

    %=======================================================
    % orthogonalize and rescale SCPMs based on changes above
    %=======================================================
    
    sol.probe.phi = orthog_modes_eigendecomp( sol.probe.phi );
    
    sol.probe.phi = enforce_scpm_fro2TOT_photonocc( sol.probe.phi, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ );

    %====================
    % Fixed probe support 
    %====================
    
%     sol.probe.support = make_rectangle( sol.sz.sz, [ 0.9 * sol.sz.r, 0.9 * sol.sz.c ]);
    sol.probe.support = make_2Dellipsoid( sol.sz.sz, [ 0.7 * sol.sz.r, 0.7 * sol.sz.c ]);

    %========
    
%     sol.probe.support = lpf_gauss( sol.probe.support, 0.10 * sol.sz.sz );
    
    %========
    
    sol.probe.support( abs( sol.probe.support ) < 1e-3 ) = 0;
    
    %========

%     sol.probe.support = ones( sol.sz.sz, 'single' );
    
    %================================================
    % load previous sample/probe/spos here if desired
    %================================================

%     tmp99 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/probes/L0105_to_L0113_combined_512x512_12Jan2022_t130245_it15001.mat' );
% 
% %     sol.probe.phi = tmp99.probe.phi;
%     
%     sol.probe.phi = [];
%     
%     sol.probe.phi( :, :, 4 ) = tmp99.probe.phi( :, :, 8 );
%     sol.probe.phi( :, :, 3 ) = tmp99.probe.phi( :, :, 7 );
%     sol.probe.phi( :, :, 2 ) = tmp99.probe.phi( :, :, 6 );
%     sol.probe.phi( :, :, 1 ) = tmp99.probe.phi( :, :, 5 );
%     
%     sol.probe.phi( :, :, 2 ) = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, 0.04 * sol.sz.sz );
%     sol.probe.phi( :, :, 1 ) = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, 0.02 * sol.sz.sz );
% 
%     sol.probe.scpm.N = size( sol.probe.phi, 3 );
   
    %========

% %     tmp0 = nocircshift2D( tmp99.probe.phi( :, :, 2 ), [-50, 20] ); 

% %     sol.probe.phi( :, :, 3 ) = tmp0;
% %     sol.probe.phi( :, :, 2 ) = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, 0.02 * sol.sz.sz );
% %     sol.probe.phi( :, :, 1 ) = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, 0.01 * sol.sz.sz );
% % 
%     sol.probe.phi = imresize3( sol.probe.phi, [ sz( 1 ), sol.sz.c, sz( 3 ) ] );
%     
% %     sol.probe.phi = sol.probe.phi( :, ( 0.5 * sz( 2 ) - 0.5 * sol.sz.c + 1 ) : ( 0.5 * sz( 2 ) + 0.5 * sol.sz.c ), : );
    

    %========================================
    % Tweak spatially incoherent probe modes?
    %========================================
    
    if ~isfield( sol.init, 'reinit_scpm_phi' ), sol.init.reinit_scpm_phi = logical( 0 ); end

    %========

    if isempty( sol.init.reinit_scpm_phi ) || sol.init.reinit_scpm_phi  

%         phi_old  = sol.probe.phi;
        sol.probe.phi = [];

    %     sol.probe.phi( :, :, 5 ) = phi_old( :, :, 3 );
    %     sol.probe.phi( :, :, 4 ) = phi_old( :, :, 2 );

    %     sol.probe.phi( :, :, 5 ) = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, 0.06 * sol.sz.sz );
    %     sol.probe.phi( :, :, 4 ) = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, 0.06 * sol.sz.sz );
        sol.probe.phi( :, :, 3 ) = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, 0.02 * sol.sz.sz );
        sol.probe.phi( :, :, 2 ) = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, 0.04 * sol.sz.sz );
        sol.probe.phi( :, :, 1 ) = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, 0.05 * sol.sz.sz );

        sol.probe.scpm.N = size( sol.probe.phi, 3 );

     %     sol.probe.phi = enforce_scpm_fro2TOT_photonocc( sol.probe.phi, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ );
         sol.probe.phi = orthog_modes_eigendecomp( sol.probe.phi );

    end

    %===========================
    % convert doubles to singles
    %===========================
    
    sol.probe.scpm.N       = single( sol.probe.scpm.N );
    sol.probe.scpm.fro2TOT = single( sol.probe.scpm.fro2TOT );
    sol.probe.scpm.max     = single( sol.probe.scpm.max );
    sol.probe.scpm.occ     = single( sol.probe.scpm.occ );
    
    %=========================
    % Shrinkwrap probe support 
    %=========================
    
    sol.swparams.blurx     = single( 0.01 ); 
    sol.swparams.blury     = single( 0.01 );
    sol.swparams.sparselvl = single( 0.60 ); 
                        
    %=======================================================
    % Sparsity constraints for probe phase ramp minimization
    %=======================================================

    sol.probe.sparse.s2DFDxy = setup_edgedetect_forwarddiff2Dxy( [], sol.sz.sz );
    
    sol.probe.sparse.support = ones( sol.sz.sz, 'logical' );  % only enforce sparsity on subregion of probe defined by some support

    sol.probe.sparse.pct_x      = single( 0.25 );
    sol.probe.sparse.pct_y      = single( 0.25 );
    sol.probe.sparse.lvl_x      = round( sol.sz.rc * sol.probe.sparse.pct_x );
    sol.probe.sparse.lvl_y      = round( sol.sz.rc * sol.probe.sparse.pct_y );
    sol.probe.sparse.threshtype = 's';
    sol.probe.sparse.threshname = 'custom_phase_ramp_min';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAMPLE INITIALIZATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %================================
    % Fixed sample support/shrinkwrap
    %================================  
 
%     sol.sample.support = 0 + 1 * make_rectangle( sol.sample.sz.sz, [ 0.43 * sol.sample.sz.r, 0.43 * sol.sample.sz.c ]);
%     % sol.probe.support = 0 + 1 * make_ellipsoid( sol.sample.sz.sz, [ 0.75 * sol.sample.sz.c, 0.75 * sol.sample.sz.r ]);
% 
%     [ sol.sample.support ] = 1 + 0 * lpf_gauss( sol.sample.support, 222.01 * sol.sample.sz.sz );
 
    %{
    
    figure; imagesc( sol.sample.support .* angle( sol.sample.T ))
    figure; imagesc( sol.sample.support )
    
    %}

    %========================================================
    % Sample magnitude and phase inequality constraint bounds
    %========================================================
 
    sol.sample.phsL = single( -0.5 * pi  );
    sol.sample.phsH = single( +0.5 * pi );

    sol.sample.absL = single( 0.00 );
    sol.sample.absH = single( 1.10 );
    
    %======================================
    % reset the sample to some random start
    %======================================
    
    if isempty( sol.sample.T ) || sol.init.reinit_sampleTF              
    
%         sol.sample.T = rand( sol.sample.sz.sz, 'single' )  + 1i * rand( sol.sample.sz.sz, 'single' );
%         sol.sample.T = lpf_gauss( sol.sample.T, 0.01 * sol.sample.sz.sz );
        
        sol.sample.T = rand( sol.sample.sz.sz, 'single' ) .* exp( 1i * 2 * pi * rand( sol.sample.sz.sz, 'single' ));
%         sol.sample.T = lpf_gauss( sol.sample.T, 0.01 * sol.sample.sz.sz );

%         sol.sample.T = sol.sample.T / max( abs( sol.sample.T( : )));
    
    end

    %=====================================================
    % Misc modifications to sample magnitude/phase scaling 
    %=====================================================

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
    
    %================================
    % Sparsity constraints for sample 
    %================================

    sol.sample.sparse.s2DFDxy = setup_edgedetect_forwarddiff2Dxy( [], sol.sample.sz.sz );
    
    sol.sample.sparse.type = 'sparse2DFDxy';
    
    sol.sample.sparse.type_sparse2DFDxy = { 'aniso_abs', ...                  % 1
                                            'aniso_phs', ...                  % 2
                                            'iso_abs', ...                    % 3
                                            'iso_phs', ...                    % 4
                                            'iso_abs_iso_phs', ...            % 5
                                            'aniso_abs_aniso_phs', ...        % 6
                                            'aniso_re_aniso_im', ...          % 7
                                            'iso_re_iso_im' };                % 8
   
    % only use sample sparsity on subregion of sample image
    sol.sample.sparse.support = ones( sol.sample.sz.sz, 'logical' ); 

    sol.sample.sparse.pct         = single( 0.25 );
    sol.sample.sparse.lvl         = round( sol.sample.sz.rc * sol.sample.sparse.pct );
    sol.sample.sparse.threshtype  = 's';
    sol.sample.sparse.threshname  = sol.sample.sparse.type_sparse2DFDxy{ 1 };
    
end

