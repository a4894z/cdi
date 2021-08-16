%

%{


%!!!!!!!!!!!!!!!!!!!!!!!!! change expt.paths.rdatah5  !!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!! change bash script         !!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!! change spos index in setup !!!!!!!!!!!!!!!!!!!!!!!!!



codelocation =  '~/Documents/Science/Matlab/Code/cdi/'; cd( codelocation ); restoredefaultpath; addpath( genpath( pwd ));



clear; close all; 

restoredefaultpath
codelocation =  '~/Documents/Science/Matlab/Code/cdi/';
cd( codelocation );
addpath( genpath( pwd ));   
clearvars -except expt sol

% clear; close all; sim_ptycho2DTPA_run

%}

%==================================================================================================

clear; close all;
restoredefaultpath; 
expt.paths.code = [ pwd, '/' ];
addpath( genpath( expt.paths.code ));
clearvars -except expt sol

rng( 'shuffle' )
reset( gpuDevice( 1 ));

%==================================================================================================

[ expt.paths ] = sim_ptycho2DTPA_paths( expt.paths ); 

% /home/ash/Documents/Science/Matlab/Code/tomolamino/testingsim_multisliceptycho.mat

%==================================================================================================

% % Create processed data:
% [ sol, expt ] = sim_ptycho2DTPA_setup( expt ); 
% return

%==================

% Load processed data
[ sol, expt ] = sim_ptycho2DTPA_loaddata( expt.paths.rdata, expt.paths );    

%==================================================================================================
% Initialize default parameters for use in phase retrieval

% sol.sample.TF   = expt.sample.TF;
% 
% sol.probe.P     = expt.probe.P;
% sol.probe.scpm  = expt.probe.scpm;

%==================

[ sol, expt ] = sim_ptycho2DTPA_defaults( sol, expt );  
  
% sol.sparse.type_sparse2DFDxy = { 'aniso_abs', ...                  % 1
%                                  'aniso_phs', ...                  % 2
%                                  'iso_abs', ...                    % 3
%                                  'iso_phs', ...                    % 4
%                                  'iso_abs_iso_phs', ...            % 5
%                                  'aniso_abs_aniso_phs', ...        % 6
%                                  'aniso_re_aniso_im', ...          % 7
%                                  'iso_re_iso_im' };                % 8     

%==================================================================================================

for ii = 1 : 5

    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    sol.phi_update = 'ER';
%     sol.sample_update = 'ePIE';
    sol.probe_update = 'DM';    % ePIE, DM
    
    sol.it.sample_sparsity = 1;
    sol.sparse.pct         = single( 0.20 );
    sol.sparse.lvl         = round( sol.sample.sz.rc * sol.sparse.pct );
    sol.sparse.threshtype  = 's';
    sol.sparse.threshname  = sol.sparse.type_sparse2DFDxy{ 7 };
    
    sol.it.spos_start = 2000e9;

    %==============
   
%     [ sol, expt ] = ptycho2DTPA_runiterations( sol, expt, 100 );      % Perform ptycho phase retrieval iterations
%     [ sol, expt ] = ptycho2DTPA_runGPUiterations_v2( sol, expt, 10 );      % Perform ptycho phase retrieval iterations
    [ sol, expt ] = ptycho2DTPA_runGPUiterations_v3( sol, expt, 10 );      % Perform ptycho phase retrieval iterations
    
    %==============
    
%     [ sol, expt ] = ptycho2DTPA_cleanup( sol, expt );            % Clean up
%     ptycho2DTPA_saveresults( expt.paths.sdata, sol, expt );      % Save results 
        
    %==============
    
    reset( gpuDevice( 1 ));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 Utility functions for run file                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ sol, expt ] = sim_ptycho2DTPA_defaults( sol, expt )

 
    warning('off','MATLAB:prnRenderer:opengl');
    
    rng( 'shuffle' )
         
    %==============================================================================================
    %========================= When to start sample/probe/spos updates ============================
    %==============================================================================================
    
    sol.it.spos_start  = 1e99;
    sol.it.spos_update = 1e99;
    
    sol.spos.suofai  = logical( 0 );                % same update order for all iterations
    sol.spos.rand_ss = 1/5;                         % 

    sol.it.sample_start        = 0;
    sol.it.sample_update       = 1;
    sol.it.sample_sparsity     = 1;
    sol.it.sample_mag_phs_ineq = 1; 
    
    sol.it.probe_start   = 25;
    sol.it.probe_update  = 1;
    sol.it.probe_scaling = 1;
    sol.it.probe_orthog  = 1;
    sol.it.probe_maxvals = 1;
    sol.it.probe_support = 1;

    sol.it.print_img_results = 50;
    sol.it.collect_metrics   = 50;
    
    %============================================
    
    sol.RAAR.beta = 0.5;
        
    %==============================================================================================
    %================================== Probe initializations =====================================
    %==============================================================================================
    
    %==================================
    % Probe initializations: Reset scpm
    %==================================
    
%     sol.probe.P( :, :, 1 ) = 0.01 + 0 * sol.probe.P( :, :, 1 );
% %     sol.probe.P( :, :, 2 ) = 0.02 + 0 * sol.probe.P( :, :, 2 );
% %     sol.probe.P( :, :, 3 ) = 0.02 + 0 * sol.probe.P( :, :, 2 );
% 
%     sol.probe.scpm.occ = [ 0.05, 0.10, 0.85 ];
%     sol.probe.scpm.occ     = sort( sol.probe.scpm.occ / norm( sol.probe.scpm.occ, 1 ));     % make sure the mode occupancy adds up to 1.0
%     
%     sol.probe.P            = enforce_scpm_fro2TOT_photonocc( sol.probe.P, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ );
%     
%     sol.probe.P            = orthog_modes_eigendecomp( sol.probe.P );

    %=================================================
    % Probe initializations: Max abs probe mode values
    %=================================================
    
    sol.probe.scpm.max = 2 * [ 70, 70, 70 ];
    
    %============================================
    % Probe initializations: Probe mode occupancy
    %============================================
    
    sol.probe.scpm.occ = [ ];
%     sol.probe.scpm.occ = [ 0.10, 0.20, 0.70 ]; 
%     sol.probe.scpm.occ = [ 0.10, 0.25, 0.65 ]; 
%     sol.probe.scpm.occ = [ 0.05, 0.15, 0.80 ];
%     sol.probe.scpm.occ = [ 0.05, 0.10, 0.85 ];
%     sol.probe.scpm.occ = [ 0.03, 0.07, 0.90 ];
    
    sol.probe.scpm.occ = sort( sol.probe.scpm.occ / norm( sol.probe.scpm.occ, 1 ));     % make sure the mode occupancy adds up to 1.0
      
    %=========================
    % probe mode total scaling
    %=========================
    
    sol.probe.scpm.fro2TOT = [];

    sol.probe.scpm.fro2TOT = 2.00 * mean( expt.meas.SI_sumD2 );
%     sol.probe.scpm.fro2TOT = 2e12;
%     sol.probe.scpm.fro2TOT = expt.probe.scpm.fro2TOT;

    %====================
    % Fixed probe support 
    %====================
    
    sol.probe.support = 0 + 1 * make_rectangle( sol.sz.sz, [ 0.7 * sol.sz.r, 0.7 * sol.sz.c ]);
%     sol.probe.support = 0 + 1 * make_2Dellipsoid( sol.sz.sz, [ 0.80 * sol.sz.r, 0.80 * sol.sz.c ]);
    [ sol.probe.support ] = 0 + 1 * lpf_gauss( sol.probe.support, 90.03 * sol.sz.sz );
    
    sol.probe.support( abs( sol.probe.support ) < 1e-3 ) = 0;
%     
    %==========================
    % Center the probe function 
    %==========================
    
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

    %==============================================================================================
    %================================= Sample initialization ======================================
    %==============================================================================================

    %================================
    % Fixed sample support/shrinkwrap
    %================================
    
%     sol.sample.support = 0 + 1 * make_rectangle( sol.sample.sz.sz, [ 0.43 * sol.sample.sz.r, 0.43 * sol.sample.sz.c ]);
%     % sol.probe.support = 0 + 1 * make_ellipsoid( sol.sample.sz.sz, [ 0.75 * sol.sample.sz.c, 0.75 * sol.sample.sz.r ]);
% 
%     [ sol.sample.support ] = 1 + 0 * lpf_gauss( sol.sample.support, 222.01 * sol.sample.sz.sz );
 
    %{
    
    figure; imagesc( sol.sample.support .* angle( sol.sample.TF ))
    figure; imagesc( sol.sample.support )
    
    %}

    %===============
    % Sample scaling 
    %===============
    
    sol.sample.phsL = -0.5 * pi;
    sol.sample.phsH = +0.5 * pi;

    sol.sample.absL = 0.00;
    sol.sample.absH = 1.00;
   
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
    
    %====================================================
    % Sample magnitude/phase scaling + misc modifications
    %====================================================
    
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

    %==============================================================================================
    %==================================== SPOS initialization =====================================
    %==============================================================================================
    
    sol.spos.update.dpx        = 1/4;  % for computing finite differences wrt scan position
    sol.spos.update.maxpx      = 4 * sol.spos.update.dpx;
    sol.spos.update.linesearch = transpose( 0.0 : sol.spos.update.dpx : sol.spos.update.maxpx );
    
%     sol.spos.update.shifttype = 'px';
    sol.spos.update.shifttype = 'subpx';

%     sol.spos.update.grad = 'all';
    sol.spos.update.grad = 'indiv';

    % for this data set, we use 0.5 um steps...only allow us to go so far from initial spos we start at
    sol.spos.update.maxcorrectr = 1.00e-6 / expt.csys.z2.dLy;
    sol.spos.update.maxcorrectc = 1.00e-6 / expt.csys.z2.dLx;
    
    sol.spos.indxsubset = expt.spos.indxsubset;                                                 % the scan positions we update the exitwaves/sample/probes over
    sol.spos.shifttype = 'px';                                                                  % When extracting part of the sample in the current scan position, use this type of pixel shifting  
                   
    %==============================================================================================
    %=========================== Gaussian LPF for reconstruction resolution =======================
    %==============================================================================================
        
    mu    = 0.5 * sol.sz.sz + 1;
    stdev = 0.80 * sol.sz.sz;
    tmp0  = make_2Dgaussian( sol.sz.sz, mu, stdev );
    
    sol.measLPF = 0 + 1 * fftshift( tmp0 );

    %==============================================================================================
    %================================ Compute initial exit waves ==================================
    %==============================================================================================
    
% %     sol.spos.N_batch  = 80;
% %     sol.spos.rs_batch = sol.spos.rs( 1 : sol.spos.N_batch, : );
% %     
% %     sol.spos.frameindx = get_indices_2Dframes( sol.spos.rs_batch, sol.sample.sz.sz, sol.sample.vs );
% % 
% %     TF = sol.sample.TF( : );
% %     TF = reshape( TF( sol.spos.frameindx ), [ sol.sz.sz, 1, sol.spos.N_batch ]);
% % 
% %     sol.phi = TF .* sol.probe.P;
%    
%     %==============
% 
%     sol.spos.frameindx = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs );
% 
%     TF = sol.sample.TF( : );
%     TF = reshape( TF( sol.spos.frameindx ), [ sol.sz.sz, 1, sol.spos.N ]);
% 
%     sol.phi = TF .* sol.probe.P;

end

%==================================================================================================
%==================================================================================================
%==================================================================================================

function [ paths ] = sim_ptycho2DTPA_paths( paths )

    %==================

    % paths.prev_sample = [ paths.code, '/G0024_apslogo_07May2020_t170828_it2201.mat' ];

    %==================

    % paths.prev_probe = [ paths.code, '/results/G0024_apslogo_27May2020_t132403_it25001.mat' ];
  
    %==================

    % paths.prev_spos = [ paths.code, '/results/G0024_apslogo_07May2020_t170828_it2201.mat' ];

    %==================

%     % path for reading & saving processed and ready-for-phase-retieval data 
%     paths.dataname = 'sim_ptycho2DTPA';
%     paths.rsdata = [ paths.code, paths.dataname, '.mat' ];

    

        % path for reading processed and ready-for-phase-retieval data 
%     paths.rdataname = '/sim_ptycho2DTPA';
    paths.rdataname = '/sim_ptycho2DTPA';
    paths.rdata = [ paths.code, paths.rdataname, '.mat' ];

    % path for saving processed and ready-for-phase-retieval data 
%     paths.sdataname = '/sim_ptycho2DTPA';
    paths.sdataname = '/sim_ptycho2DTPA';
    paths.sdata = [ paths.code, paths.sdataname, '.mat' ];

    % % path for reading & saving processed and ready-for-phase-retieval data 
    % paths.dataname = '/L0070_aps_logo';
    % paths.rsdata = [ paths.code, paths.sdataname, '.mat' ];
    
    

end
    
%==================================================================================================
%==================================================================================================
%==================================================================================================

function [ sol, expt ] = sim_ptycho2DTPA_loaddata( rdata, paths0 )

    fprintf('\n========================================================================================================'); 
    fprintf('\nLOADING Phase Retreival Experiment, \n2D Transmission Geometry and Projection Approx Assumed ...'); 

    % tmp0  = load( rsdata, 'sol' );
    % sol = tmp0.sol;
    % 
    % tmp0 = load( rsdata, 'expt' );
    % expt = tmp0.expt;

    tmp0 = load( rdata );
    sol = tmp0.sol;
    expt = tmp0.expt;
    
    expt.paths = paths0;
    
    fprintf(' done loading data!\n'); 
    fprintf('========================================================================================================\n\n'); 

    
    
%     expt.paths.code     = ldparams.codelocation;
%     expt.paths.dataname = ldparams.dataname;
%     expt.paths.rsdata   = ldparams.rsdata;

    % % check if results directory exists; if not make it
    % sol.paths.imresults = [ codelocation, 'imresults' ];
    % if exist( sol.paths.imresults, 'dir' ) ~= 7, mkdir( sol.paths.imresults ), end
    
end


%==================================================================================================
%==================================================================================================
%==================================================================================================
