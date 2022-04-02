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
    
    clear; close all; [ sol, expt] = runsolver_ptycho2DTPA;
    
    %========
    
    /net/s8iddata/export/8-id-ECA/Analysis/reduced_data/2021-3/zjiang20211214
    /net/s8iddata/export/8-id-i/2020-2/zjiang202007/reduced_data

    %}

    %================================================================================================================================================
    
    %#ok<*LOGL>
    warning('off','MATLAB:prnRenderer:opengl');
    
    rng( 'shuffle' )
    
    restoredefaultpath; 
    addpath( genpath( pwd ));
%     addpath( genpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/' ));
%     addpath( genpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/' ));
    clearvars -except expt sol
    
    %================================================================================================================================================

%     data_path = './L0105_to_L0113_combined_512x512.mat';
%     data_path = './L0258_to_L0265_combined_512x512.mat';
%     data_path = './L0266_to_L0273_combined_512x512.mat';
%     data_path = './L0274_to_L0280_combined_512x512.mat';
%     data_path = './L0302_to_L0308_combined_512x512.mat';
%     data_path = './L0333_to_L0337_combined_512x512.mat';

    %=========

% %     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/noise/sim_ptycho2DTPA.mat';      % WITH NOISE, WITHOUT BG RM

    data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/noise_rmbg/sim_ptycho2DTPA.mat'; % WITH NOISE, WITH BG RM
%     data_path = './sim_ptycho2DTPA.mat';  

    %=========
    
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p1degree_512x512.mat';  
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p11degree_512x512.mat';  
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p12degree_512x512.mat';  
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p13degree_512x512.mat';  
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p14degree_512x512.mat'; 
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p15degree_512x512.mat';  
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p16degree_512x512.mat';  
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p17degree_512x512.mat';  
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p18degree_512x512.mat'; 
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p19degree_512x512.mat'; 

%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p2degree_512x512.mat'; 
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p3degree_512x512.mat'; 
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p4degree_512x512.mat'; 
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p5degree_512x512.mat'; 
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p6degree_512x512.mat'; 
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p7degree_512x512.mat'; 
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p8degree_512x512.mat'; 
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p9degree_512x512.mat'; 
    
  

%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PetraIII_Eiffel_1p1degree_512x512.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PetraIII_Eiffel_1p11degree_512x512.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PetraIII_Eiffel_1p2degree_512x512.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PetraIII_Eiffel_1p3degree_512x512.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PetraIII_Eiffel_1p4degree_512x512.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PetraIII_Eiffel_1p5degree_512x512.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PetraIII_Eiffel_1p6degree_512x512.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PetraIII_Eiffel_1p7degree_512x512.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PetraIII_Eiffel_1p8degree_512x512.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PetraIII_Eiffel_1p9degree_512x512.mat';

%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PetraIII_Eiffel_2p0degree_256x256.mat';
    
    %=========

%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA.mat';   % WITHOUT NOISE

%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise_probe_scaling/1e-1_wrt_orig/sim_ptycho2DTPA.mat';   % WITHOUT NOISE
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise_probe_scaling/1e+1_wrt_orig/sim_ptycho2DTPA.mat';   % WITHOUT NOISE
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise_probe_scaling/1e+2_wrt_orig/sim_ptycho2DTPA.mat';   % WITHOUT NOISE
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise_probe_scaling/1e+3_wrt_orig/sim_ptycho2DTPA.mat';   % WITHOUT NOISE
    
    %================================================================================================================================================
    
%     data_path = [ pwd, '/sim_ptycho2DTPA.mat' ];
%     data_path = './pecoCSSIsimulatedmultislicedata_1024x1024.mat';

    %================================================================================================================================================
    
    load( data_path, 'sol', 'expt' );  

%     expt.paths.rsdata = data_path;
    expt.paths.rsdata = [ './', expt.paths.data_mat_name, '.mat' ];
%     expt.paths.rsdata = './pecoCSSIsimulatedmultislicedata_512x512.mat';
    
    clearvars -except expt sol
  
    %=========================================================
    % Initialize default parameters for use in phase retrieval
    %=========================================================
    
    sol.init.init_scpm_max = logical( 1 );
    sol.init.init_scpm_occ = logical( 1 );
    
    [ sol, expt ] = runsolver_ptycho2DTPA_config( sol, expt );    
    
    %================================================================================================================================================
    
    sol.rPIE_alpha_phi = single( 1e-0 );   % relatively insensitive to this?
    sol.rPIE_alpha_T   = single( 1e-1 );   
    
    %=========
    
%     sol.sample.T = rand( sol.sample.sz.sz, 'single' ) + 1i * rand( sol.sample.sz.sz, 'single' );
%     sol.sample.T = sol.sample.T / max( abs( sol.sample.T( : )));

%     % reset sample to freespace
%     sol.sample.T = ones( sol.sample, 'single' );        
         
    %=========
    % GPU Init
    %=========
    
    sol.use_gpu = true; 

    sol.gpu_id = 4; 
 
    if sol.use_gpu == true, reset( gpuDevice( sol.gpu_id )); end

    %=====================================================
    % Poisson or Gaussian Noise Model for Exitwave Updates
    %=====================================================
    
    sol.exwv_noisemodel = 'poisson';
%     sol.exwv_noisemodel = 'gaussian';
    
    sol.it.exwv_poisson_exactsearch = 1;
    
    sol.it.exwv_noisemodel_switch2poisson  = 1e99;
    sol.it.exwv_noisemodel_switch2gaussian = 1e99;
    
    %================================
    % Stochastic minibatch parameters
    %================================

    sol.spos.rand_spos_subset_pct = 0.10;      

    sol.spos.rand_spos_subset_pct = single( sol.spos.rand_spos_subset_pct );    
    
    if sol.spos.rand_spos_subset_pct > 1.00, sol.spos.rand_spos_subset_pct = single( 1.00 ); end
    
    %====================================
    % Sparse sample edges / TV parameters
    %====================================
    
%     sol.sparse.type_sparse2DFDxy = { 'aniso_abs', ...                  % 1
%                                      'aniso_phs', ...                  % 2
%                                      'iso_abs', ...                    % 3
%                                      'iso_phs', ...                    % 4
%                                      'iso_abs_iso_phs', ...            % 5
%                                      'aniso_abs_aniso_phs', ...        % 6
%                                      'aniso_re_aniso_im', ...          % 7
%                                      'iso_re_iso_im' };                % 8

    sol.sparse.pct         = single( 0.20 );
    sol.sparse.lvl         = round( sol.sample.sz.rc * sol.sparse.pct );
    sol.sparse.threshtype  = 's';
    sol.sparse.threshname  = sol.sparse.type_sparse2DFDxy{ 3 };
                              
    %======================================================
    % epoch, iteration, and iteration repeat specifications
    %======================================================
    
    N_pauseandsave = single( 2 );
    N_epochs       = single( 2500 );
    
    sol.it.probe_start     = single( 1 );      % 5 for rPIE paper
    sol.it.probe_update    = single( 1 );      % 1 for rPIE paper
    sol.it.probe_orthog    = single( 50 );     % 25 for rPIE paper
    sol.it.probe_centering = single( 1 );
    
    sol.it.sample_sparsity = single( 5e99 );
    sol.it.sample_mag_ineq = single( 1 ); 
    sol.it.sample_phs_ineq = single( 1e99 ); 
    
    sol.it.collect_metrics   = single( 10 );     % when we want to collect performance metrics
    sol.it.mkimg_meas_metric = single( 500 ); 
    sol.it.mkimg_sample_SCPM = single( 500 ); 
    
    %============
    % CHEAT CODES
    %============

%     sol.sample.T = expt.sample.T;
    
%     sol.probe.phi      = expt.probe.phi;
%     sol.probe.scpm.occ = expt.probe.scpm.occ;
%     sol.probe.scpm.N   = expt.probe.scpm.N;
%     sol.probe.scpm.max = max( reshape( abs( expt.probe.phi ), [ expt.sz.rc, expt.probe.scpm.N ] ), [], 1 );
    
%     sol.spos.rs    = expt.spos.rs;
%     sol.spos.indx  = expt.spos.indx;
    
    %================================================================================================================================================

    for ii = 1 : N_pauseandsave

        %=============
        % minibatch GD
        %=============

        [ sol, expt ] = ptycho2DTPA_runGPU_stochminibatchgrad( sol, expt, N_epochs );    

        %==============================
        % fullbatch block stochastic GD
        %==============================

%         [ sol, expt ] = ptycho2DTPA_runGPU_stochcoordgrad( sol, expt, N_epochs );  

        %===================
        % fullbatch total GD
        %===================

%         [ sol, expt ] = ptycho2DTPA_runGPU_totalgrad( sol, expt, N_epochs );     

        %==========================
        % bookkeeping and logistics
        %==========================
        
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

    %========================================
    % When to start sample/probe/spos updates 
    %========================================
    
    sol.it.spos_start  = single( 1e99 );
    sol.it.spos_update = single( 1e99 );

    sol.it.sample_start    = single( 0 );
    sol.it.sample_update   = single( 1 );
    sol.it.sample_mag_ineq = single( 1 ); 
    sol.it.sample_phs_ineq = single( 1 ); 
    sol.it.sample_sparsity = single( 1 );
    
    sol.it.probe_start     = single( 1 );      
    sol.it.probe_update    = single( 1 );
    sol.it.probe_scaling   = single( 1 );
    sol.it.probe_support   = single( 1 );
    sol.it.probe_centering = single( 1 );
    sol.it.probe_orthog    = single( 1 );
    sol.it.probe_maxvals   = single( 1 );

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

    %==================================
    % projection algorithm parameter(s)
    %==================================
    
%     sol.RAAR_beta  = single( 0.1 );
%     sol.rPIE_alpha_T = single( 1e-3 );

    %===============================================================
    % Gaussian LPF for suppressing noise at high spatial frequencies
    %===============================================================
        
    stdev = single( 0.5 * sol.sz.sz );
    tmp0  = make_2Dgaussian( sol.sz.sz, single( 0.50 * sol.sz.sz + 1 ), stdev );
    
    sol.measLPF = single( 1 + 0 * fftshift( tmp0 ));
    
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

% %     phi_old  = sol.probe.phi;
%     sol.probe.phi = [];
%     
% %     sol.probe.phi( :, :, 5 ) = phi_old( :, :, 3 );
% %     sol.probe.phi( :, :, 4 ) = phi_old( :, :, 2 );
% 
% %     sol.probe.phi( :, :, 5 ) = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, 0.06 * sol.sz.sz );
%     sol.probe.phi( :, :, 4 ) = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, 0.06 * sol.sz.sz );
%     sol.probe.phi( :, :, 3 ) = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, 0.05 * sol.sz.sz );
%     sol.probe.phi( :, :, 2 ) = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, 0.04 * sol.sz.sz );
%     sol.probe.phi( :, :, 1 ) = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, 0.03 * sol.sz.sz );
%     
%     sol.probe.scpm.N = size( sol.probe.phi, 3 );
%     
% %     sol.probe.scpm.occ = [ 0.02, 0.03, 0.95 ];
% %     sol.probe.scpm.occ = sort( sol.probe.scpm.occ / norm( sol.probe.scpm.occ, 1 ));     % make sure the mode occupancy adds up to 1.0
%      
%  %     sol.probe.phi = enforce_scpm_fro2TOT_photonocc( sol.probe.phi, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ );
%      sol.probe.phi = orthog_modes_eigendecomp( sol.probe.phi );

    %==========================
    % Max abs probe mode values
    %==========================
    
    if ~isfield( sol.probe.scpm, 'max' ), sol.probe.scpm.max = single( [ ] ); end
    
    %========
    
    if isempty( sol.probe.scpm.max ) || sol.init.init_scpm_max
    
        sol.probe.scpm.max = single( exp( +0.1 * ( 1 : sol.probe.scpm.N )));
        
        sol.probe.scpm.max = single( 200 * sol.probe.scpm.max / max( sol.probe.scpm.max( : )));
    
    end
    
    try sol.probe.scpm.max = single( reshape( sol.probe.scpm.max, [ 1, 1, sol.probe.scpm.N ] )); catch, end
    
    %=====================
    % Probe mode occupancy
    %=====================
    
    if ~isfield( sol.probe.scpm, 'occ' ), sol.probe.scpm.occ = single( [ ] ); end
    
    %========
    
    if isempty( sol.probe.scpm.occ ) || sol.init.init_scpm_occ 
        
        sol.probe.scpm.occ = transpose( exp( +1 * 0.2 * ( 1 : sol.probe.scpm.N )));                 % decaying exponential occupancy guess
        
        sol.probe.scpm.occ = single( sort( sol.probe.scpm.occ / norm( sol.probe.scpm.occ, 1 )));     % make sure the mode occupancy adds up to 1.0
        
    end

    %=========================
    % probe mode total scaling
    %=========================
    
    sol.probe.scpm.fro2TOT = single( [] );
    
%     sol.probe.scpm.fro2TOT = single( 5.0 * mean( squeeze( sum( sum( expt.meas.D .^ 2, 1 ), 2 ))));
    sol.probe.scpm.fro2TOT = single( 1.75 * mean( squeeze( sum( sum( expt.meas.D .^ 2, 1 ), 2 ))));   % for rPIE paper
%     sol.probe.scpm.fro2TOT = single( 0.9 * expt.probe.scpm.fro2TOT );

    %================================
    % orthogonalize and rescale SCPMs
    %================================
    
    sol.probe.phi = orthog_modes_eigendecomp( sol.probe.phi );
    
    sol.probe.phi = enforce_scpm_fro2TOT_photonocc( sol.probe.phi, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ );

    %====================
    % Fixed probe support 
    %====================
    
%     sol.probe.support = make_rectangle( sol.sz.sz, [ 0.9 * sol.sz.r, 0.9 * sol.sz.c ]);
    sol.probe.support = make_2Dellipsoid( sol.sz.sz, [ 0.8 * sol.sz.r, 0.8 * sol.sz.c ]);
    
%     sol.probe.support = lpf_gauss( sol.probe.support, 90.03 * sol.sz.sz );
    
    sol.probe.support( abs( sol.probe.support ) < 1e-3 ) = 0;
    
    sol.probe.support = 0 + 1 * sol.probe.support;
    
%     sol.probe.phi = sol.probe.phi .* sol.probe.support;
    
    %========
    
    sol.probe.scpm.N       = single( sol.probe.scpm.N );
    sol.probe.scpm.fro2TOT = single( sol.probe.scpm.fro2TOT );
    sol.probe.scpm.max     = single( sol.probe.scpm.max );
    sol.probe.scpm.occ     = single( sol.probe.scpm.occ );
    
    %=========================
    % Shrinkwrap probe support 
    %=========================
    
    sol.swparams.blurx     = single( 0.01 ); 
    sol.swparams.blury     = single( 0.01 );
    sol.swparams.sparselvl = single( 0.50 ); 
                        
    %===========================
    % Center the probe intensity
    %===========================
    
    %{

    figure; imagesc(abs(sol.probe.phi( :, :, end )))

    [ com ] = centerOfMass( abs( sol.probe.phi( :, :, end )));

    for pp = 1 : sol.probe.scpm.N

        sol.probe.phi( :, :, pp ) = circshift( sol.probe.phi( :, :, pp ), -1 * round( com - 0.5 * sol.sz.sz - 1 ));

    end

    figure; imagesc(abs(sol.probe.phi( :, :, end )))


    [~, Ic ] = max( sum( abs(sol.probe.phi( :, :, end )), 1 ));
    [~, Ir ] = max( sum( abs(sol.probe.phi( :, :, end )), 2 ));

    sol.probe.phi( :, :, pp ) = nocircshift2D( sol.probe.phi( :, :, pp ), -1 * round( [ Ir, Ic ] - 0.5 * sol.sz.sz - 1 ));

    figure; imagesc(abs(sol.probe.phi( :, :, end )))

    %}

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

    %===============
    % Sample scaling 
    %===============
 
    sol.sample.phsL = single( -0.5 * pi  );
    sol.sample.phsH = single( +0.5 * pi );

    sol.sample.absL = single( 0.00 );
    sol.sample.absH = single( 1.10 );

    %================================
    % Sparsity constraints for sample 
    %================================

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
    sol.sparse.threshname  = sol.sparse.type_sparse2DFDxy{ 3 };

    %====================================================
    % Sample magnitude/phase scaling + misc modifications
    %====================================================
    
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

end

