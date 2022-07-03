function [ sol, expt ] = runsolver_ptycho2DTPA( data_path )

    %#ok<*LOGL>
    warning('off','MATLAB:prnRenderer:opengl');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %{
    
    
    
    
    
    dX = convn(object,[-1,1], 'same'); 
    dY = convn(object,[-1,1]', 'same'); 
    
    

    cd /net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi
    
    restoredefaultpath; 
    clear RESTOREDEFAULTPATH_EXECUTED;
    addpath( genpath( pwd ));   
    
    %========

    clear; close all; 


    
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

%     clear; close all; % clc;

    rng( 'shuffle' )
    % rng( 555 )
    
    %==================================================================
    % Set paths to code location and add other relevant folders to path
    %==================================================================

    restoredefaultpath; 
    clear RESTOREDEFAULTPATH_EXECUTED;
    
    addpath( genpath( pwd ));
%     addpath( genpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/' ));
%     addpath( genpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/' ));
%     addpath( genpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/' ));
    addpath( genpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/' ));

    %=========================================
    % Set paths to load data from, and load it
    %=========================================
    
    if ~nargin, data_path = []; end

     % GPU2, tanzanite
    data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0341_to_L0343_combined_768x1280.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0569_to_L0570_combined_768x1280.mat';

     % GPU3, axinite
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0397_to_L0399_combined_768x1280.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0504_to_L0505_combined_768x1280.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0255_to_L0258_combined_768x1280.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0500_to_L0501_combined_768x1280.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0336_to_L0337_combined_768x1280.mat';

     % GPU3, axinite
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0496_to_L0497_combined_768x1280.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0332_to_L0333_combined_768x1280.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0492_to_L0493_combined_768x1280.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0488_to_L0489_combined_768x1280.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0324_to_L0325_combined_768x1280.mat';
 
    % GPU2, axinite
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0484_to_L0485_combined_768x1280.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0320_to_L0321_combined_768x1280.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0480_to_L0481_combined_768x1280.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0316_to_L0317_combined_768x1280.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0312_to_L0313_combined_768x1280.mat';
    
    % GPU1, axinite
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0472_to_L0473_combined_768x1280.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0308_to_L0309_combined_768x1280.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0468_to_L0469_combined_768x1280.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0304_to_L0305_combined_768x1280.mat';
%     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202206/ready_for_phaseretrieval/L0274_to_L0276_combined_768x1280.mat';

    [ sol, expt ] = runsolver_ptycho2DTPA_setpaths_loaddata( data_path );

    %==================================================================
    % Initialize/reset sample, SCPM, scan positions for phase retrieval
    %==================================================================
    
    sol.init.reinit          = logical( 1 );
    
    sol.init.reinit_misc     = logical( 1 ) & sol.init.reinit;     
    sol.init.reinit_spos     = logical( 1 ) & sol.init.reinit;
    sol.init.reinit_SCPM     = logical( 1 ) & sol.init.reinit;    
    sol.init.reinit_sampleTF = logical( 1 ) & sol.init.reinit;   

    %========

    sol.init.use_prev_spos       = logical( 0 );  
    sol.init.use_prev_SCPM       = logical( 1 );
    sol.init.use_prev_2DsampleTF = logical( 0 );

    %========
    
    if sol.init.reinit_misc,     [ sol, expt ] = define_initial_misc(          sol, expt ); end
    if sol.init.reinit_spos,     [ sol, expt ] = define_initial_scanpositions( sol, expt ); end
    if sol.init.reinit_SCPM,     [ sol, expt ] = define_initial_SCPM(          sol, expt ); end
    if sol.init.reinit_sampleTF, [ sol, expt ] = define_initial_2DsampleTF(    sol, expt ); end

    %======================================
    % !!!! CHEAT CODES FOR SIMULATIONS !!!!
    %======================================

    %     sol.sample.T = single( expt.sample.T );

%     sol.probe.phi          = single( expt.probe.phi );
%     sol.probe.scpm.occ     = single( expt.probe.scpm.occ );
%     sol.probe.scpm.N       = single( expt.probe.scpm.N );
%     sol.probe.scpm.fro2TOT = single( expt.probe.scpm.fro2TOT );
%     sol.probe.scpm.max     = single( max( reshape( abs( expt.probe.phi ), [ expt.sz.rc, expt.probe.scpm.N ] ), [], 1 ) );
% 
%     sol.spos.rs         = single( expt.spos.rs );
%     sol.spos.indxsubset = single( expt.spos.indxsubset );

    %=====================================================
    % Noise model for enforcing the measurement constraint
    %=====================================================
    
    sol.exwv_noisemodel = 'gaussian_magnitude';
%     sol.exwv_noisemodel = 'poisson';

    %=========
    % GPU Init
    %=========

    sol.use_gpu = true; 

    sol.gpu_id = 2; 

    if sol.use_gpu == true, reset( gpuDevice( sol.gpu_id )); end

    %=================================================================
    % Stochastic minibatch size ( percentage of total scan positions )
    %=================================================================

    sol.spos.rand_spos_subset_pct = 0.05;      

    C = ( sol.spos.rand_spos_subset_pct > 1.00 ) || ( sol.spos.rand_spos_subset_pct <= 0 );
    
    if C, error('Invalid Minibatch Percentage'); end

    %==================================================
    % Set default parameters for use in phase retrieval
    %==================================================
    
    sol.probe.scpm.use_fixed_occ          = logical( 0 );   % use fixed SCPM occupancies as a constraint
    sol.probe.scpm.use_fixed_scpm_scaling = logical( 1 );   % use fixed total probe intensity as a constraint
    
    sol.probe.scpm.reinit_occ          = logical( 0 );
    sol.probe.scpm.reinit_scpm_scaling = logical( 1 );
    
    [ sol, expt ] = runsolver_ptycho2DTPA_params_misc(          sol, expt );      
    [ sol, expt ] = runsolver_ptycho2DTPA_params_scanpositions( sol, expt );
    [ sol, expt ] = runsolver_ptycho2DTPA_params_SCPM(          sol, expt );
    [ sol, expt ] = runsolver_ptycho2DTPA_params_2DsampleTF(    sol, expt );

    %====================================================
    % !!!!!! MISC SCREWING AROUND WITH LOADED DATA !!!!!!
    %====================================================
    
    % testing_minibatch_spos_correction
    % return
    
% tmp0 = sol.spos.rs - min( sol.spos.rs, [], 1 );
% tmp0 = tmp0 - max( tmp0, [], 1 ) * 0.5;  
% 
% shift_px = mean( sol.spos.rs - tmp0 );
% 
% sol.spos.rs = sol.spos.rs - shift_px;
% 
% 
% sol.sample.T = nocircshift2D( sol.sample.T, shift_px );
% sol.GPU.TFvec = sol.GPU.TFvec( : );
% 
% sol.GPU.TFvec( sol.GPU.TFvec == 0 ) = 1;
%  
% for pp = 1 : sol.GPU.Nscpm
% 
%     sol.GPU.phi( :, :, pp ) = nocircshift2D( sol.GPU.phi( :, :, pp ), shift_px );
% 
% end
%         

    %===================
    % Ready, set, launch
    %===================
    
    N_pauseandsave = 2;
    N_epochs       = 250;

    for ii = 1 : N_pauseandsave

        %=============
        % Minibatch GD
        %=============

        [ sol, expt ] = ptycho2DTPA_runGPU_stochminibatchgrad( sol, expt, N_epochs );    

        %==============================
        % Fullbatch block stochastic GD
        %==============================

    %         [ sol, expt ] = ptycho2DTPA_runGPU_stochcoordgrad( sol, expt, N_epochs );    ( !!!!!!!!!! DON'T USE, NEED TO UPDATE !!!!!!!!!! )

        %===================
        % Fullbatch total GD 
        %===================

    %         [ sol, expt ] = ptycho2DTPA_runGPU_totalgrad( sol, expt, N_epochs );       ( !!!!!!!!!! DON'T USE, NEED TO UPDATE !!!!!!!!!! )

        %==========================
        % Bookkeeping and logistics
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

function [ sol, expt ] = runsolver_ptycho2DTPA_setpaths_loaddata( data_path )

    if isempty( data_path )

        %======================
        % Peco CSSI simulations
        %======================

    %     [ data_path ] = pecoCSSI_runsolver_ptycho2DTPA_setpaths_loaddata( );

        % data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p1degree_1088x2560.mat';  
        % data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p2degree_1088x2560.mat'; 
        % 
        % data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p3degree_1088x2560.mat';  
        % data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p4degree_1088x2560.mat'; 
        % 
        % data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p5degree_1088x2560.mat';  
        % data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p6degree_1088x2560.mat';  

        %=============
        % zjiang202112
        %=============

        % data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/processed_expt/L0274_to_L0280_combined_256x1024.mat';
        % data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/processed_expt/L0274_to_L0280_combined_512x1024.mat';

        % data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/processed_expt/old/L0274_to_L0280_combined_512x512.mat';
        % data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/processed_expt/old/L0105_to_L0113_combined_512x512.mat';

        %=============
        % zjiang202204
        %=============

    %     [ data_path ] = zjiang202204_runsolver_ptycho2DTPA_setpaths_loaddata( );

        %=============
        % zjiang202206
        %=============

        [ data_path ] = zjiang202206_runsolver_ptycho2DTPA_setpaths_loaddata;

        %===============
        % simulated data
        %===============

    %     data_path = [ pwd, '/sim_ptycho2DTPA.mat' ];


    %     [ data_path ] = rPIE_vs_MB_mat_runsolver_ptycho2DTPA_setpaths_loaddata( );

    %     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/noise/sim_ptycho2DTPA.mat';      % WITH NOISE, WITHOUT BG RM
    %     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/noise_rmbg/sim_ptycho2DTPA.mat'; % WITH NOISE, WITH BG RM

    %     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/noise_rmbg/sim_ptycho2DTPA_sposcorr_test.mat';  % WITH NOISE, WITH BG RM



    %     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA_sposcorr_test.mat'; % NO NOISE
    %     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA_sposcorr_shearx_plus0p5_scaley_1p2.mat'; % NO NOISE

    %     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA_sposcorr_test_rot_plus20.mat';
    %     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA_sposcorr_shearx_plus0p5_scaley_1p2.mat';
    %     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA_sposcorr_scalex1p2_scaley0p9.mat';
    %     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA_sposcorr_shearx_minux0p2_sheary_plus0p4.mat';



    %     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA_sposcorr_test.mat'; % NO NOISE
    %     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA_sposcorr_shearx_plus0p5_scaley_1p2.mat'; % NO NOISE

    %     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA_sposcorr_test_rot_plus20.mat';
    %     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA_sposcorr_shearx_plus0p5_scaley_1p2.mat';
    %     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA_sposcorr_scalex1p2_scaley0p9.mat';
    %     data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA_sposcorr_shearx_minux0p2_sheary_plus0p4.mat';



    %     data_path = './sim_ptycho2DTPA_sposcorr_test_plus20rot.mat';
    %     data_path = './sim_ptycho2DTPA_sposcorr_shearx_minux0p2_sheary_plus0p4.mat';
    %     data_path = './sim_ptycho2DTPA_sposcorr_scalex1p2_scaley0p9.mat';
    %     data_path = './sim_ptycho2DTPA_sposcorr_shearx_plus0p5_scaley_1p2.mat';

    %     data_path = './sim_ptycho2DTPA_sposcorr_test.mat';


        % data_path = './L0105_to_L0113_combined_512x512.mat';
    %     data_path = './L0078_to_L0081_combined_512x512.mat';
    %     data_path = './L0339_to_L0341_combined_512x512.mat';

        % data_path = './L0274_to_L0280_combined_512x512.mat';
        % data_path = './L0020_to_L0032_combined_512x512.mat';

    %     data_path = './L0342_to_L0345_combined_1024x512.mat';

    end

    %========================
    % Load the processed data 
    %========================

    % load( data_path, 'expt' );  
    load( data_path, 'sol', 'expt' );  

    %====================================
    % Update path for currently read data
    %====================================

    %     expt.paths.rsdata = data_path;
    expt.paths.rsdata = [ './', expt.paths.data_mat_name, '.mat' ];
    %     expt.paths.rsdata = './pecoCSSIsimulatedmultislicedata.mat';

    clearvars -except expt sol
    
    if ~exist( 'sol', 'var' ), sol = struct; end
    
end

%====================================================================================================================================================

function [ sol, expt ] = runsolver_ptycho2DTPA_params_misc( sol, expt )     

    %================================================================
    % Poisson noise model step length parameters for exitwave updates
    %================================================================

    if strcmp( sol.exwv_noisemodel, 'poisson' )
    
%         sol.poissonexwv.steplength_update_type = 'use_exact_linesearch';
%         sol.poissonexwv.steplength_update_type = 'use_sign_test';   % we should always initially use an exact line search? and we should just choose from these two?
        
        sol.poissonexwv.steplength_update_type = 'use_recursive';
    %     sol.poissonexwv.steplength_update_type = 'use_previous';
    
        %========
        
        sol.poissonexwv.steplength_update_freq = 1;   
        
        %========
                
        sol.poissonexwv.alpha_start          = single( 0.50 );
        sol.poissonexwv.delta_alpha_signtest = single( 0.05 );

        %========

        sol.poissonexwv.Nalpha = single( 5 );

        %========
        
%         sol.poissonexwv.alpha_max = exp( -0.2 * ( 1 : sol.probe.scpm.N ));
%         sol.poissonexwv.alpha_max = sol.poissonexwv.alpha_max - min( sol.poissonexwv.alpha_max );
%         sol.poissonexwv.alpha_max = sol.poissonexwv.alpha_max / max( sol.poissonexwv.alpha_max );
%         sol.poissonexwv.alpha_max = single( 3 * sol.poissonexwv.alpha_max + 1 );

    %     sol.poissonexwv.alpha_max = 6 - 0.65 * ( 1 : sol.probe.scpm.N );
    %     sol.poissonexwv.alpha_max( sol.poissonexwv.alpha_max < 0 ) = 0.5;

        sol.poissonexwv.alpha_max            = 2 * ones( 1, sol.probe.scpm.N, 'single' );
        sol.poissonexwv.alpha_max( end )     = 1;
        sol.poissonexwv.alpha_max( end - 1 ) = 1;
        
        %========
        
        sol.poissonexwv.alpha_minmax         = zeros( 2, sol.probe.scpm.N, 'single' );
        sol.poissonexwv.alpha_minmax( 2, : ) = sol.poissonexwv.alpha_max;

    %     figure; plot( sol.poissonexwv.alpha_max )

    end
    
    %=====================================
    % rPIE adaptive step length parameters
    %=====================================
    
    sol.rPIE_alpha_phi = single( 1e-0 );   % relatively insensitive to this?
    sol.rPIE_alpha_T   = single( 2.5e-1 );   

    %========================================
    % When to start sample/probe/spos updates 
    %========================================
    
    % sol.it --> sol.constraints
    
    sol.it.spos_start  = single( 250 + 1 );
    sol.it.spos_update = single( 10 );

    sol.it.sample_start     = single( 0 );
    sol.it.sample_update    = single( 1 );
    sol.it.sample_mag_ineq  = single( 1 ); 
    sol.it.sample_phs_ineq  = single( 1 ); 
    sol.it.sample_sparsity  = single( 1 );
%     sol.it.sample_smoothing = single( 1 );
    
    sol.it.probe_start        = single( 250 + 1 );      
    sol.it.probe_update       = single( 10 );
    sol.it.probe_scaling      = single( 1 );
    sol.it.probe_support      = single( 1 );
    sol.it.probe_centering    = single( 1 );
    sol.it.probe_orthog       = single( 25 );   % 25 for rPIE paper
    sol.it.probe_maxvals      = single( 1 );
    sol.it.probe_smoothing    = single( 1 );
    sol.it.probe_phase_sparse = single( 1e99 );
    
    sol.it.collect_metrics   = single( 25 );     % when we want to collect performance metrics
    sol.it.mkimg_meas_metric = single( 100 ); 
    sol.it.mkimg_sample_SCPM = single( 50 ); 

    %========

    % for metrics plotting
    
    sol.metrics.metrics_ylim      = [ 10^-4, 10^0 ];
    sol.metrics.grad_metrics_ylim = [ 10^-5, 10^0 ];
    
    sol.metrics.legend_loc = 'northeast';
%     sol.metrics.legend_loc = 'southwest';

    %===============================================================
    % Gaussian LPF for suppressing noise at high spatial frequencies
    %===============================================================
        
    sol.measLPF = make_2Dgaussian( sol.sz.sz, single( 0.50 * sol.sz.sz + 1 ), single( 0.25 * sol.sz.sz ) );
    
    sol.measLPF = single( fftshift( sol.measLPF / max( abs( sol.measLPF( : )))));

    
end

%====================================================================================================================================================

function [ sol, expt ] = runsolver_ptycho2DTPA_params_SCPM( sol, expt )   
    
    %========================
    % Guess of SCPM occupancy
    %========================
   
    if sol.probe.scpm.use_fixed_occ || sol.probe.scpm.reinit_occ
        
        sol.probe.scpm.occ = transpose( exp( +1 * 1.2 * ( 1 : sol.probe.scpm.N )));         % decaying exponential occupancy guess
        sol.probe.scpm.occ = sort( sol.probe.scpm.occ / norm( sol.probe.scpm.occ, 1 ));     % make sure the mode occupancy adds up to 1.0
        
    else
        
        sol.probe.scpm.occ = [];
        
    end
    
    %=======================================================
    % Guess of probe photons (sum of probe intensity pixels)
    %=======================================================

    if sol.probe.scpm.use_fixed_scpm_scaling || sol.probe.scpm.reinit_scpm_scaling
    
        sol.probe.scpm.fro2TOT = single( 5.0 * mean( squeeze( sum( sum( expt.meas.D .^ 2, 1 ), 2 ))));    
%         sol.probe.scpm.fro2TOT = single( ( 1 + 0.05 * ( 2 * rand - 1 )) * expt.probe.scpm.fro2TOT ); 
%         sol.probe.scpm.fro2TOT = single( 1.02 * expt.probe.scpm.fro2TOT ); 
        
    else
        
        sol.probe.scpm.fro2TOT = [];
        
    end
        
    %====================
    % Orthogonalize SCPMs
    %====================

%     [ sol.probe.phi ] = orthog_modes_eigendecomp( sol.probe.phi );

    %==============
    % Rescale SCPMs
    %==============

    [ sol.probe.phi ] = enforce_scpm_fro2TOT_photonocc( sol.probe.phi, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ );
    
    %================================
    % Max abs probe mode pixel values 
    %================================

    % ?????????????? IS THIS NEEDED? CAN THIS BE DONE AWAY WITH BY SCPM SMOOTHING ?????????????? 
    % ?????????????? OR GAUSSIAN METRIC W/ LPF ?????????????? 

    sol.probe.scpm.max = [];

    sol.probe.scpm.max = single( exp( +0.25 * ( 1 : sol.probe.scpm.N )));
    sol.probe.scpm.max = single( 2250 * sol.probe.scpm.max / max( sol.probe.scpm.max( : )));

%     sol.probe.scpm.max = single( [ 30, 60, 160 ]);

    %============================
    % Gaussian smoothing of SCPMs
    %============================
    
    sol.probe.gaussian_blur_std_x = 0.10 * sol.sz.c;
    sol.probe.gaussian_blur_std_y = 0.10 * sol.sz.r;
    
    sol.probe.gaussian_lpf = make_2Dgaussian( sol.sz.sz, 0.5 * sol.sz.sz + 1, [ sol.probe.gaussian_blur_std_y, sol.probe.gaussian_blur_std_x ] );
    
    sol.probe.gaussian_lpf = fftshift( sol.probe.gaussian_lpf );
    
    %=========================
    % Shrinkwrap probe support 
    %=========================
    
    sol.probe.swparams_blur_x    = single( 0.01 ); 
    sol.probe.swparams_blur_y    = single( 0.01 );
    sol.probe.swparams_sparselvl = single( 0.60 ); 

    %============================
    % Fixed probe support options
    %============================
    
    sol.probe.support = make_rectangle( sol.sz.sz, [ 0.90 * sol.sz.r, 0.90 * sol.sz.c ]);
%     sol.probe.support = make_2Dellipsoid( sol.sz.sz, [ 0.80 * sol.sz.r, 0.80 * sol.sz.c ]);
    
%     sol.probe.support = lpf_gauss( sol.probe.support, 0.01 * sol.sz.sz );

    sol.probe.support( abs( sol.probe.support ) < 1e-3 ) = 0;

    %=======================================================
    % Sparsity constraints for probe phase ramp minimization
    %=======================================================

    sol.probe.sparse.s2DFDxy = setup_edgedetect_forwarddiff2Dxy( [], sol.sz.sz );
    
    sol.probe.sparse.support = ones( sol.sz.sz, 'logical' );  % only enforce sparsity on subregion of probe defined by some support

    sol.probe.sparse.pct_x      = single( 0.50 );
    sol.probe.sparse.pct_y      = single( 0.50 );
    sol.probe.sparse.lvl_x      = round( sol.sz.rc * sol.probe.sparse.pct_x );
    sol.probe.sparse.lvl_y      = round( sol.sz.rc * sol.probe.sparse.pct_y );
    sol.probe.sparse.threshtype = 's';
    sol.probe.sparse.threshname = 'xy_aniso_phase';

    %===========================
    % convert doubles to singles
    %===========================

    sol.probe.phi          = single( sol.probe.phi );
    sol.probe.scpm.N       = single( sol.probe.scpm.N );
    sol.probe.scpm.fro2TOT = single( sol.probe.scpm.fro2TOT );
    sol.probe.scpm.occ     = single( sol.probe.scpm.occ );
    sol.probe.scpm.max     = single( sol.probe.scpm.max );
    
end

%====================================================================================================================================================

function [ sol, expt ] = runsolver_ptycho2DTPA_params_scanpositions( sol, expt )   

%     sol.spos.indxsubset = single( sol.spos.indx );   % FOR OLD SIMULATED DATA
    
    sol.spos.rs         = single( sol.spos.rs );
    sol.spos.indxsubset = single( sol.spos.indxsubset );
        
    %========================================================
    % keep track of the initial scan positions used in ptycho
    %========================================================
    
    if ~isfield( sol.spos, 'rs0' ), sol.spos.rs0 = expt.spos.rs( sol.spos.indxsubset, : ); end
    
    %========================================================================================
    % Gradient descent based scan position correction w/ scan positions as decision variables
    %========================================================================================

    sol.spos.correction.Naalpha_rs = 4;  
    sol.spos.correction.aalpha_rs  = linspace( 0.0, 3.00, sol.spos.correction.Naalpha_rs );   

    sol.spos.correction.optimize_rs_GD = 'individual';
%     sol.spos.correction.optimize_rs_GD = 'collective';

    sol.spos.correction.noise_model = sol.exwv_noisemodel;
    % sol.spos.correction.noise_model = 'poisson';
%     sol.spos.correction.noise_model = 'gaussian_magnitude';

    sol.spos.correction.maxcorrect_r = single( 250 );
    sol.spos.correction.maxcorrect_c = single( 250 );
    
    sol.spos.correction.reset_rand_r = single( 5 );
    sol.spos.correction.reset_rand_c = single( 5 );

    sol.spos.correction.rs0          = sol.spos.rs0;

    %=================================================================================================
    % When extracting part of the sample in the current scan position, use this type of pixel shifting  
    %=================================================================================================
    
    sol.spos.shifttype = 'px';   
    
    %============================================================================
    % center the scan positions around the origin, shift sample/probe accordingly
    %============================================================================
    
    %{
    
    tmp0 = sol.spos.rs - min( sol.spos.rs, [], 1 );
    tmp0 = tmp0  - max( tmp0 , [], 1 ) * 0.5;  

    shift_yx = sum( tmp0 - sol.spos.rs ) / sol.spos.N;

    sol.spos.rs = sol.spos.rs + shift_yx;
    
    sol.sample.T = nocircshift2D( sol.sample.T, -1 * round( double( shift_yx ))) ;
    
%     for pp = 1 : sol.probe.scpm.N,  sol.probe.phi( :, :, pp ) = nocircshift2D( sol.probe.phi( :, :, pp ), -shift_yx ); end
    
    %}
    
    %============================================
    % Modifications to the current scan positions
    %============================================
    
    % RANDOM
    
%     sol.spos.rs = sol.spos.rs + 10 * random( sol.spos.N, 2, 'single' );

    %========
    
%     % ROTATION
%     
%     theta = +20 - 10;
%     sol.spos.rotation = [ [ +cosd( theta ), +sind( theta )  ]; ...
%                           [ -sind( theta ), +cosd( theta )  ] ];
%     
%     sol.spos.rs = transpose( sol.spos.rotation * transpose( sol.spos.rs ));

    %========
        
%     % SHEAR
%     
%     s_x = +0.5;
%     sol.spos.shear_x = [ [ 1,   0 ]; ...
%                          [ s_x, 1 ] ];
%     
%     sol.spos.rs = transpose( sol.spos.shear_x * transpose( sol.spos.rs ));
% 
%     s_y = +0.00;
%     sol.spos.shear_y = [ [ 1, s_y ]; ...
%                          [ 0, 1   ] ];
%                      
%     sol.spos.rs = transpose( sol.spos.shear_y * transpose( sol.spos.rs ));
    
    %========
    
%     % SCALING
%     
%     sol.spos.rs( :, 1 ) = 1.1 * sol.spos.rs( :, 1 );   % scaling in rows ( y )
%     sol.spos.rs( :, 2 ) = 1.0 * sol.spos.rs( :, 2 );   % scaling in cols ( x )

    



end

%====================================================================================================================================================

function [ sol, expt ] = runsolver_ptycho2DTPA_params_2DsampleTF( sol, expt )   

    %================================
    % Fixed sample support/shrinkwrap
    %================================  
 
%     sol.sample.support = 0 + 1 * make_rectangle( sol.sample.sz.sz, [ 0.43 * sol.sample.sz.r, 0.43 * sol.sample.sz.c ]);
%     % sol.probe.support = 0 + 1 * make_ellipsoid( sol.sample.sz.sz, [ 0.75 * sol.sample.sz.c, 0.75 * sol.sample.sz.r ]);
% 
%     [ sol.sample.support ] = 1 + 0 * lpf_gauss( sol.sample.support, 222.01 * sol.sample.sz.sz );
%  
    %{
    
    figure; imagesc( sol.sample.support .* angle( sol.sample.T ))
    figure; imagesc( sol.sample.support )
    
    %}

    %========================================================
    % Sample magnitude and phase inequality constraint bounds
    %========================================================
 
    sol.sample.phsL = single( -0.9 * pi  );
    sol.sample.phsH = single( +0.9 * pi );

    sol.sample.absL = single( 0.00 );
    sol.sample.absH = single( 1.00 );
  
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

    sol.sample.sparse.pct         = single( 0.40 );
    sol.sample.sparse.lvl         = round( sol.sample.sz.rc * sol.sample.sparse.pct );
    sol.sample.sparse.threshtype  = 's';
    sol.sample.sparse.threshname  = sol.sample.sparse.type_sparse2DFDxy{ 1 };
    

end

%====================================================================================================================================================


