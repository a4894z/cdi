%
%{


%!!!!!!!!!!!!!!!!!!!!!!!!! change expt.paths.rdatah5  !!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!! change bash script         !!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!! change spos index in setup !!!!!!!!!!!!!!!!!!!!!!!!!



codelocation = '~/Documents/Science/Matlab/Code/cdi/'; cd( codelocation ); restoredefaultpath; addpath( genpath( pwd ));



clear; close all; 

restoredefaultpath
codelocation =  '~/Documents/Science/Matlab/Code/cdi/';
cd( codelocation );
addpath( genpath( pwd ));   
clearvars -except expt sol

% clear; close all; zjiang202011_L0002_apslogo_200nm_run

%}

%==================================================================================================

clear; close all;
% restoredefaultpath; 
expt.paths.code = pwd;
addpath( genpath( expt.paths.code ));
clearvars -except expt sol

rng( 'shuffle' )
reset( gpuDevice( 1 ));

% try
%     
%    gpuArray(1);
%    canUseGPU=true;
%    
% catch
%     
%  canUseGPU=false;
%  
% end

%==================================================================================================

[ expt.paths ] = zjiang202011_L0002_apslogo_200nm_paths( expt.paths ); 

%==================================================================================================

process_ptycho_data = logical( 0 ); %#ok<LOGL>

if process_ptycho_data == true
    
%     expt.paths.rdatah5 = '/home/8ididata/2020-3/zjiang202011/reduced_data/Avg_L0002_apslogo_200nm_0001.h5';

    [ sol, expt ] = zjiang202011_L0002_apslogo_200nm_setup( expt );              % ***** Create ***** processed data:
    
else
    
    
    [ sol, expt ] = zjiang202007_loaddata( expt.paths.sdata, expt.paths );       % ***** Load ***** processed data

end

%==================================================================================================
% Initialize default parameters for use in phase retrieval

[ sol, expt ] = zjiang202011_L0002_apslogo_200nm_defaults( sol, expt );
  
% sol.sparse.type_sparse2DFDxy = { 'aniso_abs', ...                  % 1
%                                  'aniso_phs', ...                  % 2
%                                  'iso_abs', ...                    % 3
%                                  'iso_phs', ...                    % 4
%                                  'iso_abs_iso_phs', ...            % 5
%                                  'aniso_abs_aniso_phs', ...        % 6
%                                  'aniso_re_aniso_im', ...          % 7
%                                  'iso_re_iso_im' };                % 8

%==================================================================================================

for ii = 1 : 100
    
    %============================================
    
%     sol.phi_update = 'RAAR';
%     sol.RAAR.beta = 1/2;
% 
% %     sol.sample_update = 'ePIE';
%     sol.probe_update = 'DM';    % ePIE, DM, mrPIE, rPIE
% 
%     sol.it.sample_sparsity = 1;
%     sol.sparse.pct         = single( 0.04 );
%     sol.sparse.lvl         = round( sol.sample.sz.rc * sol.sparse.pct );
%     sol.sparse.threshtype  = 'h';
%     sol.sparse.threshname  = sol.sparse.type_sparse2DFDxy{ 1 };
% 
%     sol.it.spos_start = 1e99;
% 
%     %==============
%     
%     [ sol, expt ] = ptycho2DTPA_runiterations( sol, expt, 200 );      % Perform ptycho phase retrieval iterations                          
%           
%     %==============
%     
%     [ sol, expt ] = ptycho2DTPA_cleanup( sol, expt );          % Clean up             
%     ptycho2DTPA_saveresults( expt.paths.sdata, sol, expt );    % Save results    
    
    %============================================
    
    sol.phi_update = 'ER';
%     sol.sample_update = 'ePIE';
    sol.probe_update = 'DM';    % ePIE, DM
    
    sol.it.sample_sparsity = 1e0;
    sol.sparse.pct         = single( 0.10 );
    sol.sparse.lvl         = round( sol.sample.sz.rc * sol.sparse.pct );
    sol.sparse.threshtype  = 's';
    sol.sparse.threshname  = sol.sparse.type_sparse2DFDxy{ 3 };
    
    sol.it.spos_start = 2000e9;

    %==============
   
%     [ sol, expt ] = ptycho2DTPA_runiterations( sol, expt, 100 );      % Perform ptycho phase retrieval iterations
%     [ sol, expt ] = ptycho2DTPA_runGPUiterations( sol, expt, 10 );      % Perform ptycho phase retrieval iterations
    [ sol, expt ] = ptycho2DTPA_runGPUiterations_v3( sol, expt, 10 );      % Perform ptycho phase retrieval iterations
    
    %==============
    
    [ sol, expt ] = ptycho2DTPA_cleanup( sol, expt );          % Clean up             
    ptycho2DTPA_saveresults( expt.paths.sdata, sol, expt );    % Save results 
    
    %============================================
    
%     sol.phi_update = 'ER';
% %     sol.sample_update = 'ePIE';
%     sol.probe_update = 'DM';    % ePIE, DM
%     
%     sol.it.sample_sparsity = 1;
%     sol.sparse.pct         = single( 0.20 );
%     sol.sparse.lvl         = round( sol.sample.sz.rc * sol.sparse.pct );
%     sol.sparse.threshtype  = 's';
%     sol.sparse.threshname  = sol.sparse.type_sparse2DFDxy{ 1 };
%     
%     sol.it.spos_start = 2000e9;
% 
%     %==============
%    
% %     [ sol, expt ] = ptycho2DTPA_runiterations( sol, expt, 100 );      % Perform ptycho phase retrieval iterations
% %     [ sol, expt ] = ptycho2DTPA_runGPUiterations( sol, expt, 10 );      % Perform ptycho phase retrieval iterations
% %     [ sol, expt ] = ptycho2DTPA_runGPUiterations_v2( sol, expt, 10 );      % Perform ptycho phase retrieval iterations
%     [ sol, expt ] = ptycho2DTPA_runGPUiterations_v3( sol, expt, 10 );      % Perform ptycho phase retrieval iterations
%     
%     %==============
%     
%     [ sol, expt ] = ptycho2DTPA_cleanup( sol, expt );          % Clean up             
%     ptycho2DTPA_saveresults( expt.paths.sdata, sol, expt );    % Save results    
        
    %============================================
    
    reset( gpuDevice( 1 ));
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------- Utility functions for run file -----------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ sol, expt ] = zjiang202011_L0002_apslogo_200nm_defaults( sol, expt )     

    warning('off','MATLAB:prnRenderer:opengl');
    
    rng( 'shuffle' )
         
    %==============================================================================================
    %------------------------- When to start sample/probe/spos updates ----------------------------
    %==============================================================================================
    
    sol.it.spos_start  = 1e99;
    sol.it.spos_update = 1e99;
    
    sol.spos.suofai  = logical( 0 );                % same update order for all iterations
    sol.spos.rand_ss = 2/5;                         % 

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

    sol.it.print_img_results = 100;
    sol.it.collect_metrics   = 100;
    
    %============================================
    
    sol.RAAR.beta = 0.5;
        
    %==============================================================================================
    %----------------------------------- Probe initializations ------------------------------------
    %==============================================================================================
    
    % Here, define probe constraint reconstruction parameters
    
    %============================================
    % Probe initializations: Reset scpm
    %============================================
    
    sol.probe.P( :, :, 1 ) = 0.01 + 0 * sol.probe.P( :, :, 1 );
    sol.probe.P( :, :, 2 ) = 0.02 + 0 * sol.probe.P( :, :, 2 );
%     sol.probe.P( :, :, 3 ) = 0.02 + 0 * sol.probe.P( :, :, 3 );

    sol.probe.scpm.occ = [ 0.05, 0.10, 0.85 ];
    sol.probe.scpm.occ = sort( sol.probe.scpm.occ / norm( sol.probe.scpm.occ, 1 ));     % make sure the mode occupancy adds up to 1.0
    
%     sol.probe.P = enforce_scpm_fro2TOT_photonocc( sol.probe.P, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ );
    sol.probe.P = orthog_modes_eigendecomp( sol.probe.P );

    %============================================
    % Probe initializations: Max abs probe mode values
    %============================================
    
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
      
%     sol.probe.P = orthog_modes_eigendecomp( sol.probe.P );

    %============================================
    % probe mode total scaling
    %============================================
    
%     sol.probe.scpm.fro2TOT = [];
    sol.probe.scpm.fro2TOT = 15.00 * mean( expt.meas.SI_sumD2 );
%     sol.probe.scpm.fro2TOT = 5e9;

%     sol.probe.P = enforce_scpm_fro2TOT_photonocc( sol.probe.P, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ );
    
    %============================================
    % Fixed probe support 
    %============================================
    
    sol.probe.support = 0 + 1 * make_rectangle( sol.sz.sz, [ 0.5 * sol.sz.r, 0.7 * sol.sz.c ]);
    % sol.probe.support = 0 + 1 * make_ellipsoid( sol.sz.sz, [ 0.75 * sol.sz.c, 0.75 * sol.sz.r ]);
    [ sol.probe.support ] = 1 + 0 * lpf_gauss( sol.probe.support, 90.03 * sol.sz.sz );
    
    sol.probe.support( abs( sol.probe.support ) < 1e-3 ) = 0;
    
%     sol.probe.P = sol.probe.P .* sol.probe.support;
    
    %============================================
    % Center the probe function 
    %============================================
    
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
    %----------------------------------- Sample initialization ------------------------------------
    %==============================================================================================
    %
    % Fixed sample support/shrinkwrap
    %
    %============================================
    
%     sol.sample.support = 0 + 1 * make_rectangle( sol.sample.sz.sz, [ 0.43 * sol.sample.sz.r, 0.43 * sol.sample.sz.c ]);
%     % sol.probe.support = 0 + 1 * make_ellipsoid( sol.sample.sz.sz, [ 0.75 * sol.sample.sz.c, 0.75 * sol.sample.sz.r ]);
% 
%     [ sol.sample.support ] = 1 + 0 * lpf_gauss( sol.sample.support, 222.01 * sol.sample.sz.sz );
 
    %{
    
    figure; imagesc( sol.sample.support .* angle( sol.sample.TF ))
    figure; imagesc( sol.sample.support )
    
    %}

    %============================================
    %
    % Sample scaling 
    %
    %============================================
    
    sol.sample.phsL = -0.20 * pi;
    sol.sample.phsH = +0.20 * pi;

    sol.sample.absL = 0.0;
    sol.sample.absH = 1.0;
   
    %============================================
    %
    % Sparsity constraints for sample 
    %
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
   
    % only use sample sparsity on subregion of sample image
    sol.sparse.support = ones( sol.sample.sz.sz, 'logical' ); 
                                 
    %============================================
    %
    % Sample magnitude/phase scaling + misc modifications
    %
    %============================================
    
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
    %------------------------------------ SPOS initialization -------------------------------------
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
                   
    %================================================================
    %--------- Gaussian LPF for reconstruction resolution -----------
    %================================================================
        
    mu    = 0.5 * sol.sz.sz + 1;
    stdev = 0.90 * sol.sz.sz;
    tmp0  = make_2Dgaussian( sol.sz.sz, mu, stdev );
    
    sol.measLPF = 0 + 1 * fftshift( tmp0 );

    %============================================
    %--------- Compute initial exit waves  ------
    %============================================
    
%     sol.spos.frameindx = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs );
% 
%     TF = sol.sample.TF( : );
%     TF = reshape( TF( sol.spos.frameindx ), [ sol.sz.sz, 1, sol.spos.N ]);
% %     TF = reshape( TF, [ sol.sz.sz, 1, sol.spos.N ] );
% 
%     sol.phi = TF .* sol.probe.P;
% 
%     clear( 'TF' )

end


%==================================================================================================

function [ paths ] = zjiang202011_L0002_apslogo_200nm_paths( paths )

    % read hdf5 data path
    paths.rdatah5 = '/media/ash/Saltwater/Data/experimental/zjiang202011/cssi/reduced_data/Avg_L0002_apslogo_200nm_0001.h5';

    %==================

    % % paths.prev_sample = [ paths.code, '' ];
%     paths.prev_sample = [ paths.code, '/results/L0367_A4rot2Dptycho000deg_12Aug2020_t145041_it3866.mat' ];
    paths.prev_sample = [ paths.code, '/results/Avg_L0002_apslogo_200nm_0001_11Dec2020_t074748_it60001.mat' ];
    
    %==================

    paths.prev_probe = [ paths.code, '/results/Avg_L0002_apslogo_200nm_0001_11Dec2020_t074748_it60001.mat' ];
%     paths.prev_probe = [ paths.code, '/results/Avg_L0002_apslogo_200nm_0001_09Dec2020_t022324_it140001.mat' ];
%     paths.prev_probe = [ paths.code, '/results/ptycho0078_apslogo_125nm_15Nov2020_t014518_it20001.mat' ];

    %==================

    % paths.prev_spos = [ paths.code, '' ];

    %==================

    % path for reading processed and ready-for-phase-retieval data 
    paths.rdataname = '/Avg_L0002_apslogo_200nm_0001';
    paths.rdata = [ paths.code, paths.rdataname, '.mat' ];

    % path for saving processed and ready-for-phase-retieval data 
    paths.sdataname = '/Avg_L0002_apslogo_200nm_0001';
    paths.sdata = [ paths.code, paths.sdataname, '.mat' ];

    %==================
    
    % the known missing data blemish due to dead pixels, detector segments, beamstop, etc 

%     paths.blemish = '~/Documents/Science/Matlab/Code/cdi/run/missingdetectorpixels/blemish_lambda750k.mat';
    paths.blemish = '/home/beams/ATRIPATH/Documents/MATLAB/Code/cdi/run/missingdetectorpixels/blemish_lambda750k.mat';
    
end

%==================================================================================================
