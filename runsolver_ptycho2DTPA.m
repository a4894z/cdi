%
%#ok<*LOGL>
warning('off','MATLAB:prnRenderer:opengl');

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


7.5 x10^9
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; % clc;

rng( 'shuffle' )
% rng( 555 )

%==============================================================
% Set paths to code location and add relevant folders and files
%==============================================================

restoredefaultpath; 
addpath( genpath( pwd ));
% addpath( genpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/' ));
% addpath( genpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/' ));
addpath( genpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/' ));
clearvars -except expt sol

%============================
% Set paths to load data from
%============================

% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/noise/sim_ptycho2DTPA.mat';      % WITH NOISE, WITHOUT BG RM
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/noise_rmbg/sim_ptycho2DTPA.mat'; % WITH NOISE, WITH BG RM

data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA_sposcorr_test.mat'; % NO NOISE

%=========
% 
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p1degree_1088x2560.mat';  
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p2degree_1088x2560.mat'; 
% 
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p3degree_1088x2560.mat';  
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p4degree_1088x2560.mat'; 
% 
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p5degree_1088x2560.mat';  
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/PETRAIII_Eiffel_Small_noTi_74ol_1p6degree_1088x2560.mat';  

%=========

% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/processed_expt/L0274_to_L0280_combined_256x1024.mat';
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/processed_expt/L0274_to_L0280_combined_512x1024.mat';

% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/processed_expt/old/L0274_to_L0280_combined_512x512.mat';
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/processed_expt/old/L0105_to_L0113_combined_512x512.mat';

%=========

% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/processed_expt/L0020_to_L0032_combined_512x512.mat';   % +0
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/processed_expt/L0020_to_L0032_combined_512x1024.mat';  % +0
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/processed_expt/L0062_to_L0065_combined_1024x128.mat';  % +60
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/processed_expt/L0058_to_L0061_combined_512x1024.mat';  % +30
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/processed_expt/L0078_to_L0081_combined_512x512.mat';   % +180
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/processed_expt/L0089_to_L0092_combined_1152x128.mat';  % +30
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/processed_expt/L0107_to_L0112_combined_768x256.mat';   % +30
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/processed_expt/L0113_to_L0116_combined_768x1024.mat';  % +5
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/processed_expt/L0117_to_L0120_combined_768x1024.mat';  % -5 
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/processed_expt/L0130_to_L0135_combined_1024x512.mat';  % -30
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/processed_expt/L0156_to_L0159_combined_512x512.mat';   % +10
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/processed_expt/L0160_to_L0163_combined_512x512.mat';   % -10

% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0145_to_L0153_combined_1024x256.mat';   % +(?) 90

% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0336_to_L0338_combined_512x512.mat';  % +0
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0339_to_L0341_combined_512x512.mat';  % +180

% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0342_to_L0345_combined_1024x512.mat';  % +5
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0346_to_L0349_combined_1024x512.mat';  % -5
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0350_to_L0353_combined_1024x512.mat';  % +10
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0354_to_L0357_combined_1024x512.mat';  % -10

% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0358_to_L0362_combined_1024x256.mat';  % +15
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0363_to_L0367_combined_1024x256.mat';  % -15 

% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0368_to_L0372_combined_1024x256.mat';  % +20
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0373_to_L0377_combined_1024x256.mat';  % -20

% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0390_to_L0395_combined_1024x256.mat';  % +30
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0396_to_L0401_combined_1024x256.mat';  % -30

% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0416_to_L0422_combined_1024x256.mat';  % +40
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0423_to_L0429_combined_1024x256.mat';  % -40

% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0430_to_L0437_combined_1024x256.mat';  % +45
% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0438_to_L0445_combined_1024x256.mat';  % -45

% data_path = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202204/ready_for_phaseretrieval/L0510_to_L0518_combined_1024x128.mat';  % +90

%=========


% data_path = './sim_ptycho2DTPA_sposcorr_test.mat';

% data_path = [ pwd, '/sim_ptycho2DTPA.mat' ];
% data_path = './L0105_to_L0113_combined_512x512.mat';

% data_path = './L0274_to_L0280_combined_512x512.mat';
% data_path = './L0020_to_L0032_combined_512x512.mat';

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

%==================================================================
% Initialize/reset sample, SCPM, scan positions for phase retrieval
%==================================================================

sol.init.reinit_spos     = logical( 0 );   
sol.init.reinit_SCPM     = logical( 0 );   
sol.init.reinit_sampleTF = logical( 0 );
sol.init.reinit_misc     = logical( 0 );   

sol.init.use_prev_spos       = logical( 0 );  
sol.init.use_prev_SCPM       = logical( 1 );
sol.init.use_prev_2DsampleTF = logical( 0 );

if sol.init.reinit_misc,     [ sol, expt ] = define_initial_misc(          sol, expt ); end
if sol.init.reinit_spos,     [ sol, expt ] = define_initial_scanpositions( sol, expt ); end
if sol.init.reinit_SCPM,     [ sol, expt ] = define_initial_SCPM(          sol, expt ); end
if sol.init.reinit_sampleTF, [ sol, expt ] = define_initial_2DsampleTF(    sol, expt ); end

%==================================================
% Set default parameters for use in phase retrieval
%==================================================

[ sol, expt ] = runsolver_ptycho2DTPA_params_misc(          sol, expt );      
[ sol, expt ] = runsolver_ptycho2DTPA_params_scanpositions( sol, expt );
[ sol, expt ] = runsolver_ptycho2DTPA_params_SCPM(          sol, expt );
[ sol, expt ] = runsolver_ptycho2DTPA_params_2DsampleTF(    sol, expt );

%=========
% GPU Init
%=========

sol.use_gpu = true; 

sol.gpu_id = 1; 

if sol.use_gpu == true, reset( gpuDevice( sol.gpu_id )); end

%================================
% Stochastic minibatch parameters
%================================

sol.spos.rand_spos_subset_pct = 0.25;      

sol.spos.rand_spos_subset_pct = single( sol.spos.rand_spos_subset_pct );    

if sol.spos.rand_spos_subset_pct > 1.00, sol.spos.rand_spos_subset_pct = single( 1.00 ); end

%======================================
% !!!! CHEAT CODES FOR SIMULATIONS !!!!
%======================================

%     sol.sample.T = expt.sample.T;

    sol.probe.phi          = expt.probe.phi;
    sol.probe.scpm.occ     = expt.probe.scpm.occ;
    sol.probe.scpm.N       = expt.probe.scpm.N;
    sol.probe.scpm.fro2TOT = expt.probe.scpm.fro2TOT;
    sol.probe.scpm.max     = max( reshape( abs( expt.probe.phi ), [ expt.sz.rc, expt.probe.scpm.N ] ), [], 1 );
    
%     sol.spos.rs    = expt.spos.rs0;
%     sol.spos.indx  = expt.spos.indx;
    
%=============================================
% Misc last minute stuff before running ptycho
%=============================================

%     sol.probe.scpm.max = single( exp( +0.1 * ( 1 : sol.probe.scpm.N )));
%     sol.probe.scpm.max = 2.5 * single( 200 * sol.probe.scpm.max / max( sol.probe.scpm.max( : )));
% 
%     sol.probe.scpm.fro2TOT = 2.5 * single( 2.00 * mean( squeeze( sum( sum( expt.meas.D .^ 2, 1 ), 2 )))); 
% 
%     sol.probe.phi = orthog_modes_eigendecomp( sol.probe.phi );
%     sol.probe.phi = enforce_scpm_fro2TOT_photonocc( sol.probe.phi, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ );
 
%========

% testing_minibatch_spos_correction
% return

%================================================================================================================================================
    
N_pauseandsave = single( 1 );
N_epochs       = single( 200 );

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

%====================================================================================================================================================

function [ sol, expt ] = runsolver_ptycho2DTPA_params_misc( sol, expt )     

    %=====================================================
    % Poisson or Gaussian Noise Model for Exitwave Updates
    %=====================================================
 
   sol.exwv_noisemodel = 'gaussian';
    
%      sol.exwv_noisemodel = 'poisson';
%      sol.it.exwv_poisson_exactsearch = 5;
    
    %===================================================================
    % hold the SCPM occupancies and total scaling fixed, or reinitialize
    %===================================================================
    
    sol.probe.scpm.use_fixed_occ = logical( 0 );   
    sol.probe.scpm.reinit_occ    = logical( 0 );
    
    sol.probe.scpm.use_fixed_scpm_scaling = logical( 0 );  
    sol.probe.scpm.reinit_scpm_scaling    = logical( 0 );

    %=====================================
    % rPIE adaptive step length parameters
    %=====================================
    
    sol.rPIE_alpha_phi = single( 1e-0 );   % relatively insensitive to this?
    sol.rPIE_alpha_T   = single( 1e-1 );   

    %========================================
    % When to start sample/probe/spos updates 
    %========================================
    
    sol.it.spos_start  = single( 5 );
    sol.it.spos_update = single( 5 );

    sol.it.sample_start     = single( 0 );
    sol.it.sample_update    = single( 1 );
    sol.it.sample_mag_ineq  = single( 1 ); 
    sol.it.sample_phs_ineq  = single( 1e99 ); 
    sol.it.sample_sparsity  = single( 1e99 );
%     sol.it.sample_smoothing = single( 1 );
    
    sol.it.probe_start        = single( 400e99 + 1 );      
    sol.it.probe_update       = single( 1 );
    sol.it.probe_scaling      = single( 1 );
    sol.it.probe_support      = single( 1 );
    sol.it.probe_centering    = single( 1 );
    sol.it.probe_orthog       = single( 50 );   % 25 for rPIE paper
    sol.it.probe_maxvals      = single( 1 );
    sol.it.probe_smoothing    = single( 1 );
    sol.it.probe_phase_sparse = single( 1e99 );
    
    sol.it.collect_metrics   = single( 20 );     % when we want to collect performance metrics
    sol.it.mkimg_meas_metric = single( 20 ); 
    sol.it.mkimg_sample_SCPM = single( 20 ); 
    
    %=============================================
    % Performance metric and epoch update counters
    %=============================================
    
    if ~isfield( sol.it, 'mtot' ),  sol.it.mtot  = 1; end
    if ~isfield( sol.it, 'metr' ),  sol.it.metr  = 1; end
    if ~isfield( sol.it, 'epoch' ), sol.it.epoch = 1; end

    %========

    % for metrics plotting
    
    sol.metrics.metrics_ylim      = [ 10^-4, 10^0 ];
    sol.metrics.grad_metrics_ylim = [ 10^-5, 10^0 ];
    
    sol.metrics.legend_loc = 'northeast';
%     sol.metrics.legend_loc = 'southwest';
    
    %========
    
    % when monitoring the Poisson cost function metric, we need a scaling offset so that it's always positive
    
%     I_m                        = expt.meas.D .^ 2;
%     log_I_m                    = log( I_m );
%     log_I_m( isinf( log_I_m )) = 0;
% 
%     sol.metrics.poisson_offset = ( I_m - I_m .* log_I_m );
%     sol.metrics.poisson_offset = sum( sol.metrics.poisson_offset(:) ) / size( I_m, 3 );
% 
%     clear( 'I_m', 'log_I_m' )

    %===============================================================
    % Gaussian LPF for suppressing noise at high spatial frequencies
    %===============================================================
        
    sol.measLPF = make_2Dgaussian( sol.sz.sz, single( 0.50 * sol.sz.sz + 1 ), single( 0.5 * 1e5 * sol.sz.sz ) );
    
    sol.measLPF = fftshift( sol.measLPF / max( abs( sol.measLPF( : ))));

%     expt.meas.D = sol.measLPF .* expt.meas.D;
    
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
    
end

%====================================================================================================================================================

function [ sol, expt ] = runsolver_ptycho2DTPA_params_SCPM( sol, expt )   

    %========================
    % Guess of SCPM occupancy
    %========================
   
    if sol.probe.scpm.use_fixed_occ || sol.probe.scpm.reinit_occ
        
        sol.probe.scpm.occ = transpose( exp( +1 * 1.8 * ( 1 : sol.probe.scpm.N )));         % decaying exponential occupancy guess
        sol.probe.scpm.occ = sort( sol.probe.scpm.occ / norm( sol.probe.scpm.occ, 1 ));     % make sure the mode occupancy adds up to 1.0
        
    else
        
        sol.probe.scpm.occ = [];
        
    end
    
    %=======================================================
    % Guess of probe photons (sum of probe intensity pixels)
    %=======================================================

    if sol.probe.scpm.use_fixed_scpm_scaling || sol.probe.scpm.reinit_scpm_scaling
    
        sol.probe.scpm.fro2TOT = single( 9.5 * mean( squeeze( sum( sum( expt.meas.D .^ 2, 1 ), 2 )))); 
        
    else
        
        sol.probe.scpm.fro2TOT = [];
        
    end
        
    %==============
    % Orthogonalize
    %==============

%     [ sol.probe.phi ] = orthog_modes_eigendecomp( sol.probe.phi );
% 
%     %========
%     % Rescale
%     %========
% 
%     [ sol.probe.phi ] = enforce_scpm_fro2TOT_photonocc( sol.probe.phi, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ );
    
    %================================
    % Max abs probe mode pixel values 
    %================================

    % ?????????????? IS THIS NEEDED? CAN THIS BE DONE AWAY WITH BY SCPM SMOOTHING ?????????????? 
    % ?????????????? OR GAUSSIAN METRIC W/ LPF ?????????????? 

    sol.probe.scpm.max = [];

    sol.probe.scpm.max = single( exp( +0.65 * ( 1 : sol.probe.scpm.N )));
    sol.probe.scpm.max = single( 120 * sol.probe.scpm.max / max( sol.probe.scpm.max( : )));

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
    sol.probe.swparams_sparselvl = single( 0.80 ); 

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
    
    sol.spos.rs         = single( sol.spos.rs );
    sol.spos.indxsubset = single( sol.spos.indxsubset );
    
    %================================================
    % Gradient descent based scan position correction
    %================================================



sol.spos.correction.Naalpha_rs = 6;  
sol.spos.correction.aalpha_rs  = linspace( 0.0, 5.00, sol.spos.correction.Naalpha_rs );   

% sol.spos.correction.optimize_rs_GD = 'individual';
sol.spos.correction.optimize_rs_GD = 'collective';

% sol.spos.correction.noise_model = 'poisson';
sol.spos.correction.noise_model = 'gaussian';
                


    %=================================================================================================
    % When extracting part of the sample in the current scan position, use this type of pixel shifting  
    %=================================================================================================
    
    sol.spos.shifttype = 'px';   
    
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
 
    sol.sample.phsL = single( -0.5 * pi  );
    sol.sample.phsH = single( +0.5 * pi );

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

    sol.sample.sparse.pct         = single( 0.5 );
    sol.sample.sparse.lvl         = round( sol.sample.sz.rc * sol.sample.sparse.pct );
    sol.sample.sparse.threshtype  = 's';
    sol.sample.sparse.threshname  = sol.sample.sparse.type_sparse2DFDxy{ 3 };
    

end

%====================================================================================================================================================


