%

%{

clear; close all; testing_mb_vs_rPIEalpha

%}

%====================================================================================================================================================

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/';

%=============================================
% define the paths for the raw data to analyze
%=============================================

[ path_data, N_trials ] = define_paths( rootpath_data );

path_data_fields = fieldnames( path_data );

%=============================================================
% load previously processed metrics variable, do some plotting
%=============================================================

% load /net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/metrics_no_noise_MB_vs_rPIEalpha.mat;
% 
% %========
% 
% plot_data_IDs = { 'mb_0p05', ...
%                   'mb_0p10', ...
%                   'mb_0p20', ...
%                   'mb_0p33', ...
%                   'mb_0p50' };
%               
% plot_data_IDs = path_data_fields;
%       
% %========
% 
% N_alphatrials   = 10;
% N_metricsepochs = 51;
% 
% metrics_plotting( length( plot_data_IDs ) ) = struct;
% 
% for jj = 1 : length( plot_data_IDs )                                    % loop over the minibatch tests
%     
%     metrics_plotting( jj ).data_id = ( plot_data_IDs{ jj } );
% 
%     for kk = 1 : length( path_data.( plot_data_IDs{ jj } ) )            % loop over the 10 trials
%         
%         metrics_plotting( jj ).rPIEalpha{ kk } = metrics.( plot_data_IDs{ jj } ){ kk }.rPIE_alpha( 1 );
%     
%         metrics_plotting( jj ).log10_meas_all_avg{ kk } = metrics.( plot_data_IDs{ jj } ){ kk }.log10_meas_all_avg;
%         metrics_plotting( jj ).log10_meas_all_max{ kk } = metrics.( plot_data_IDs{ jj } ){ kk }.log10_meas_all_max;
%         metrics_plotting( jj ).log10_meas_all_min{ kk } = metrics.( plot_data_IDs{ jj } ){ kk }.log10_meas_all_min;
%         metrics_plotting( jj ).log10_meas_all_std{ kk } = std( metrics.( plot_data_IDs{ jj } ){ kk }.log10_meas_all, [], 2 );
% 
%     end
% 
% end
% 
% 
% return

%====================================================================
% load and process the raw results, save to separate metrics variable
%====================================================================

for jj = 1 : length( path_data_fields )
 
    for kk = 1 : length( path_data.( path_data_fields{ jj } ) )
        
        metrics.( path_data_fields{ jj } ){ kk } = load_and_plot( path_data.( path_data_fields{ jj } ){ kk }, ...
                                                                  N_trials.( path_data_fields{ jj } )( kk ),  ...
                                                                  path_data_fields{ jj } );                       
    
    end

end

%====================================================================================================================================================

function metrics = load_and_plot( path_data, N_trials, path_data_fields )

    y_lim = [-1, 6];
    
    %========
    
    sim_ptycho2DTPA = cell( N_trials, 1 );
    
    for ii = 1 : N_trials

        sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
        
    end

    %========

    meas_all_avg = 0;
    
    for ii = 1 : N_trials
        
        log10_meas_all( :, ii ) = log10( sim_ptycho2DTPA{ ii }.sol.metrics.meas_all ); %#ok<AGROW>

        meas_all_avg = meas_all_avg + log10_meas_all( :, ii );

    end

    meas_all_avg = meas_all_avg / N_trials;

    log10_meas_all_min = transpose( min( log10_meas_all, [], 2 ));
    log10_meas_all_max = transpose( max( log10_meas_all, [], 2 ));  

%     name_data = num2str( [ sim_ptycho2DTPA{ 1 }.sol.rPIE_alpha, sim_ptycho2DTPA{ 1 }.sol.spos.rand_spos_subset_pct ], 'rPIE_alpha = %0.8f, MBpct = %0.4f');

    %=================================================
    % create a struct from metrics data for future use
    %=================================================
    
    for ii = 1 : N_trials
        
        metrics.meas_all( :, ii )          = sim_ptycho2DTPA{ ii }.sol.metrics.meas_all;
        metrics.timing( ii )               = sim_ptycho2DTPA{ ii }.sol.timings;
        metrics.it( ii )                   = sim_ptycho2DTPA{ ii }.sol.it;
        metrics.rPIE_alpha( ii )           = sim_ptycho2DTPA{ ii }.sol.rPIE_alpha;
        
        if isfield( sim_ptycho2DTPA{ ii }.sol.spos, 'rand_spos_subset_pct' )
            
            metrics.rand_spos_subset_pct( ii ) = sim_ptycho2DTPA{ ii }.sol.spos.rand_spos_subset_pct;
            
        end
        
    end
    
    %========
    
    if isfield( sim_ptycho2DTPA{ ii }.sol.spos, 'rand_spos_subset_pct' )
        
        name_data = [ path_data_fields, num2str( [ metrics.rPIE_alpha( 1 ), metrics.rand_spos_subset_pct( 1 ) ], ', rPIE_alpha = %0.9f, MBpct = %0.4f') ];
        save_data = [ path_data_fields, num2str( [ metrics.rPIE_alpha( 1 ), metrics.rand_spos_subset_pct( 1 ) ], '_rPIEalpha_%0.9f_MBpct_%0.4f'), '.jpg' ];
        
    else
        
        name_data = [ path_data_fields, num2str( sim_ptycho2DTPA{ii}.sol.rPIE_alpha, ', rPIE_alpha = %0.9f') ];
        save_data = [ path_data_fields, num2str( sim_ptycho2DTPA{ii}.sol.rPIE_alpha, '_rPIEalpha_%0.9f'), '.jpg' ];
        
    end
    
    %========
    
    metrics.log10_meas_all_max  = log10_meas_all_max;
    metrics.log10_meas_all_min  = log10_meas_all_min;
    metrics.log10_meas_all      = log10_meas_all;
    metrics.log10_meas_all_avg  = meas_all_avg;
    metrics.name_data           = name_data;
    metrics.N_trials            = N_trials;
    metrics.path_data           = path_data;
    
    %============================================================
    % create figures of the specified metrics over the trials run
    %============================================================
    
%     skip = 1;
%     
%     h1 = figure();  
%     set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
% 
%     x     = metrics.it( 1 ).mtot( 1 : skip : end );
%     xconf = [ x, x( end : -1 : 1) ]; 
%     yconf = [ log10_meas_all_max, fliplr( log10_meas_all_min ) ];
%     
%     p = fill( xconf, yconf, 'red' );
%     p.FaceColor = [ 1, 0.8, 0.8 ];      
%     p.EdgeColor = 'none';           
% 
%     hold on
% 
%     plot( metrics.it( 1 ).mtot( 1 : skip : end ), meas_all_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
%     xlabel('Epoch')
%     ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
%     hold off
%     title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
%     grid on
%     ylim( y_lim )
% %     legend( 'Location', 'northeast' ) 
%     
%     hold off
     
    %========
    
    skip = 1;

    h1 = figure();  
    set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )

    hold on

    for ii = 1 : N_trials

        plot( metrics.it( ii ).mtot( 1 : skip : end ), log10( metrics.meas_all( 1 : skip : end, ii )  ), '-', 'linewidth', 2 )

    end

    plot( metrics.it( 1 ).mtot( 1 : skip : end ), meas_all_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )

    xlabel('Epoch')
    ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
    hold off
    title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
    grid on
    ylim( y_lim )
    legend( 'Location', 'northeast' ) 
    
    hold off

    export_fig( save_data, '-r120.0' )
    close all;

end

%====================================================================================================================================================

function [ path_data, N_trials ] = define_paths( rootpath_data )

%================
% NO NOISE MB 10%
%================

path_data.mb_0p01 = {};
N_trials.mb_0p01  = [];

path_data.mb_0p01{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p01_oldsample_for_probe_update/independenttrials_20Sep2021_t143405/' ];
N_trials.mb_0p01( end + 1 )  = 10;

path_data.mb_0p01{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p001_oldsample_for_probe_update/independenttrials_20Sep2021_t144026/' ];
N_trials.mb_0p01( end + 1 )  = 10;

path_data.mb_0p01 = transpose( path_data.mb_0p01 );
N_trials.mb_0p01  = transpose( N_trials.mb_0p01 );

%=============
% NOISY MB 10%
%=============

path_data.mb_0p01_noise = {};
N_trials.mb_0p01_noise  = [];

path_data.mb_0p01_noise{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p01_noise/independenttrials_15Nov2021_t144212/' ];
N_trials.mb_0p01_noise( end + 1 )  = 10;

path_data.mb_0p01_noise{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p001_noise/independenttrials_15Nov2021_t230419/' ];
N_trials.mb_0p01_noise( end + 1 )  = 10;

path_data.mb_0p01_noise = transpose( path_data.mb_0p01_noise );
N_trials.mb_0p01_noise  = transpose( N_trials.mb_0p01_noise );

return

%===============
% Stoch Block GD
%===============

path_data.stoch_bgd = {};
N_trials.stoch_bgd  = [];

path_data.stoch_bgd{ end + 1 } = [ rootpath_data, 'blind_block_stoch/cdi_rPIE_stochGD_alpha1p00_oldsample_for_probe_update/independenttrials_10Oct2021_t140518/' ];
N_trials.stoch_bgd( end + 1 )  = 10;

path_data.stoch_bgd{ end + 1 } = [ rootpath_data, 'blind_block_stoch/cdi_rPIE_stochGD_alpha0p10_oldsample_for_probe_update/independenttrials_03Oct2021_t161143/' ];
N_trials.stoch_bgd( end + 1 )  = 10;

path_data.stoch_bgd{ end + 1 } = [ rootpath_data, 'blind_block_stoch/cdi_rPIE_stochGD_alpha0p01_oldsample_for_probe_update/independenttrials_22Sep2021_t104744/' ];
N_trials.stoch_bgd( end + 1 )  = 10;

path_data.stoch_bgd{ end + 1 } = [ rootpath_data, 'blind_block_stoch/cdi_rPIE_stochGD_alpha0p001_oldsample_for_probe_update/independenttrials_12Nov2021_t192759/' ];
N_trials.stoch_bgd( end + 1 )  = 10;

path_data.stoch_bgd{ end + 1 } = [ rootpath_data, 'blind_block_stoch/cdi_rPIE_stochGD_alpha0p0001_oldsample_for_probe_update/independenttrials_13Nov2021_t113330/'];
N_trials.stoch_bgd( end + 1 )  = 10;

path_data.stoch_bgd{ end + 1 } = [ rootpath_data, 'blind_block_stoch/cdi_rPIE_stochGD_alpha0p00001_oldsample_for_probe_update/independenttrials_13Nov2021_t182247/' ];
N_trials.stoch_bgd( end + 1 )  = 10;

path_data.stoch_bgd{ end + 1 } = [ rootpath_data, 'blind_block_stoch/cdi_rPIE_stochGD_alpha0p000001_oldsample_for_probe_update/independenttrials_13Nov2021_t022040/' ];
N_trials.stoch_bgd( end + 1 )  = 10;

path_data.stoch_bgd{ end + 1 } = [ rootpath_data, 'blind_block_stoch/cdi_rPIE_stochGD_alpha0p0000001_oldsample_for_probe_update/independenttrials_12Nov2021_t081029/'];
N_trials.stoch_bgd( end + 1 )  = 10;

path_data.stoch_bgd{ end + 1 } = [ rootpath_data, 'blind_block_stoch/cdi_rPIE_stochGD_alpha0p00000001_oldsample_for_probe_update/independenttrials_11Nov2021_t112402/' ];
N_trials.stoch_bgd( end + 1 )  = 10;

path_data.stoch_bgd{ end + 1 } = [ rootpath_data, 'blind_block_stoch/cdi_rPIE_stochGD_alpha0p000000001_oldsample_for_probe_update/independenttrials_10Nov2021_t122727/' ];
N_trials.stoch_bgd( end + 1 )  = 10;

path_data.stoch_bgd = transpose( path_data.stoch_bgd );
N_trials.stoch_bgd  = transpose( N_trials.stoch_bgd );

%=============
% 1% Minibatch
%=============

path_data.mb_0p01 = {};
N_trials.mb_0p01  = [];

path_data.mb_0p01{ end + 1 } = [ rootpath_data, 'mb0p01/cdi_rPIE_mb0p01_alpha1p00_oldsample_for_probe_update/independenttrials_05Nov2021_t143406/' ];
N_trials.mb_0p01( end + 1 )  = 10;

path_data.mb_0p01{ end + 1 } = [ rootpath_data, 'mb0p01/cdi_rPIE_mb0p01_alpha0p10_oldsample_for_probe_update/independenttrials_13Nov2021_t023404/' ];
N_trials.mb_0p01( end + 1 )  = 10;

path_data.mb_0p01{ end + 1 } = [ rootpath_data, 'mb0p01/cdi_rPIE_mb0p01_alpha0p01_oldsample_for_probe_update/independenttrials_01Oct2021_t161445/' ];
N_trials.mb_0p01( end + 1 )  = 10;

path_data.mb_0p01{ end + 1 } = [ rootpath_data, 'mb0p01/cdi_rPIE_mb0p01_alpha0p001_oldsample_for_probe_update/independenttrials_05Nov2021_t143333/' ];
N_trials.mb_0p01( end + 1 )  = 10;

path_data.mb_0p01{ end + 1 } = [ rootpath_data, 'mb0p01/cdi_rPIE_mb0p01_alpha0p0001_oldsample_for_probe_update/independenttrials_06Nov2021_t105931/' ];
N_trials.mb_0p01( end + 1 )  = 10;

path_data.mb_0p01{ end + 1 } = [ rootpath_data, 'mb0p01/cdi_rPIE_mb0p01_alpha0p00001_oldsample_for_probe_update/independenttrials_07Nov2021_t062512/' ];
N_trials.mb_0p01( end + 1 )  = 10;

path_data.mb_0p01{ end + 1 } = [ rootpath_data, 'mb0p01/cdi_rPIE_mb0p01_alpha0p000001_oldsample_for_probe_update/independenttrials_08Nov2021_t024900/' ];
N_trials.mb_0p01( end + 1 )  = 10;

path_data.mb_0p01{ end + 1 } = [ rootpath_data, 'mb0p01/cdi_rPIE_mb0p01_alpha0p0000001_oldsample_for_probe_update/independenttrials_08Nov2021_t230507/' ];
N_trials.mb_0p01( end + 1 )  = 10;

path_data.mb_0p01{ end + 1 } = [ rootpath_data, 'mb0p01/cdi_rPIE_mb0p01_alpha0p00000001_oldsample_for_probe_update/independenttrials_10Nov2021_t012243/' ];
N_trials.mb_0p01( end + 1 )  = 10;

path_data.mb_0p01{ end + 1 } = [ rootpath_data, 'mb0p01/cdi_rPIE_mb0p01_alpha0p000000001_oldsample_for_probe_update/independenttrials_11Nov2021_t213037/' ];
N_trials.mb_0p01( end + 1 )  = 10;

path_data.mb_0p01 = transpose( path_data.mb_0p01 );
N_trials.mb_0p01  = transpose( N_trials.mb_0p01 );

%=============
% 5% Minibatch
%=============

path_data.mb_0p05 = {};
N_trials.mb_0p05  = [];

path_data.mb_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha1p00_oldsample_for_probe_update/independenttrials_30Oct2021_t051441/' ];
N_trials.mb_0p05( end + 1 )  = 10;

path_data.mb_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p10_oldsample_for_probe_update/independenttrials_23Oct2021_t210230/' ];
N_trials.mb_0p05( end + 1 )  = 10;

path_data.mb_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p01_oldsample_for_probe_update/independenttrials_24Sep2021_t123331/' ];
N_trials.mb_0p05( end + 1 )  = 10;

path_data.mb_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p001_oldsample_for_probe_update/independenttrials_04Nov2021_t042443/' ];
N_trials.mb_0p05( end + 1 )  = 10;

path_data.mb_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p0001_oldsample_for_probe_update/independenttrials_06Nov2021_t105744/' ];
N_trials.mb_0p05( end + 1 )  = 10;

path_data.mb_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p00001_oldsample_for_probe_update/independenttrials_06Nov2021_t215004/' ];
N_trials.mb_0p05( end + 1 )  = 10;

path_data.mb_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p000001_oldsample_for_probe_update/independenttrials_07Nov2021_t073855/' ];
N_trials.mb_0p05( end + 1 )  = 10;

path_data.mb_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p0000001_oldsample_for_probe_update/independenttrials_07Nov2021_t182814/' ];
N_trials.mb_0p05( end + 1 )  = 10;

path_data.mb_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p00000001_oldsample_for_probe_update/independenttrials_08Nov2021_t052208/' ];
N_trials.mb_0p05( end + 1 )  = 10;

path_data.mb_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p000000001_oldsample_for_probe_update/independenttrials_08Nov2021_t161250/' ];
N_trials.mb_0p05( end + 1 )  = 10;

path_data.mb_0p05 = transpose( path_data.mb_0p05 );
N_trials.mb_0p05  = transpose( N_trials.mb_0p05 );

%==============
% 10% Minibatch
%==============

path_data.mb_0p10 = {};
N_trials.mb_0p10  = [];

path_data.mb_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha1p00_oldsample_for_probe_update/independenttrials_23Oct2021_t190740/' ];
N_trials.mb_0p10( end + 1 )  = 10;

path_data.mb_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p10_oldsample_for_probe_update/independenttrials_10Nov2021_t140344/'];
N_trials.mb_0p10( end + 1 )  = 10;

path_data.mb_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p01_oldsample_for_probe_update/independenttrials_20Sep2021_t143405/' ];
N_trials.mb_0p10( end + 1 )  = 10;

path_data.mb_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p001_oldsample_for_probe_update/independenttrials_20Sep2021_t144026/' ];
N_trials.mb_0p10( end + 1 )  = 10;

path_data.mb_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p0001_oldsample_for_probe_update/independenttrials_01Nov2021_t113502/' ];
N_trials.mb_0p10( end + 1 )  = 10;

path_data.mb_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p00001_oldsample_for_probe_update/independenttrials_01Nov2021_t212739/' ];
N_trials.mb_0p10( end + 1 )  = 10;

path_data.mb_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p000001_oldsample_for_probe_update/independenttrials_02Nov2021_t071445/' ];
N_trials.mb_0p10( end + 1 )  = 10;

path_data.mb_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p0000001_oldsample_for_probe_update/independenttrials_02Nov2021_t175702/' ];
N_trials.mb_0p10( end + 1 )  = 10;

path_data.mb_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p00000001_oldsample_for_probe_update/independenttrials_03Nov2021_t045630/' ];
N_trials.mb_0p10( end + 1 )  = 10;

path_data.mb_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p000000001_oldsample_for_probe_update/independenttrials_03Nov2021_t155907/' ];
N_trials.mb_0p10( end + 1 )  = 10;

path_data.mb_0p10 = transpose( path_data.mb_0p10 );
N_trials.mb_0p10  = transpose( N_trials.mb_0p10 );

%==============
% 20% Minibatch
%==============

path_data.mb_0p20 = {};
N_trials.mb_0p20  = [];

path_data.mb_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha1p00_oldsample_for_probe_update/independenttrials_19Oct2021_t191344/' ];
N_trials.mb_0p20( end + 1 )  = 10;

path_data.mb_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p10_oldsample_for_probe_update/independenttrials_12Oct2021_t073340/' ];
N_trials.mb_0p20( end + 1 )  = 10;

path_data.mb_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p01_oldsample_for_probe_update/independenttrials_21Sep2021_t205815/' ];
N_trials.mb_0p20( end + 1 )  = 10;

path_data.mb_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p001_oldsample_for_probe_update/independenttrials_21Sep2021_t205346/' ];
N_trials.mb_0p20( end + 1 )  = 10;

path_data.mb_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p0001_oldsample_for_probe_update/independenttrials_24Sep2021_t121129/' ];
N_trials.mb_0p20( end + 1 )  = 10;

path_data.mb_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p00001_oldsample_for_probe_update/independenttrials_24Sep2021_t224259/' ];
N_trials.mb_0p20( end + 1 )  = 10;

path_data.mb_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p000001_oldsample_for_probe_update/independenttrials_25Sep2021_t091618/' ];
N_trials.mb_0p20( end + 1 )  = 10;

path_data.mb_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p0000001_oldsample_for_probe_update/independenttrials_25Sep2021_t195104/' ];
N_trials.mb_0p20( end + 1 )  = 10;

path_data.mb_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p00000001_oldsample_for_probe_update/independenttrials_01Nov2021_t113429/' ];
N_trials.mb_0p20( end + 1 )  = 10;

path_data.mb_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p000000001_oldsample_for_probe_update/independenttrials_01Nov2021_t225858/' ];
N_trials.mb_0p20( end + 1 )  = 10;

path_data.mb_0p20 = transpose( path_data.mb_0p20 );
N_trials.mb_0p20  = transpose( N_trials.mb_0p20 );

%==============
% 33% Minibatch
%==============

path_data.mb_0p33 = {};
N_trials.mb_0p33  = [];

path_data.mb_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha1p00_oldsample_for_probe_update/independenttrials_07Oct2021_t204113/' ];
N_trials.mb_0p33( end + 1 )  = 10;

path_data.mb_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p10_oldsample_for_probe_update/independenttrials_03Oct2021_t155218/' ];
N_trials.mb_0p33( end + 1 )  = 10;

path_data.mb_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p01_oldsample_for_probe_update/independenttrials_25Sep2021_t132343/' ];
N_trials.mb_0p33( end + 1 )  = 10;

path_data.mb_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p001_oldsample_for_probe_update/independenttrials_26Sep2021_t002144/' ];
N_trials.mb_0p33( end + 1 )  = 10;

path_data.mb_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p0001_oldsample_for_probe_update/independenttrials_26Sep2021_t112108/' ];
N_trials.mb_0p33( end + 1 )  = 10;

path_data.mb_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p00001_oldsample_for_probe_update/independenttrials_03Nov2021_t145655/' ];
N_trials.mb_0p33( end + 1 )  = 10;

path_data.mb_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p000001_oldsample_for_probe_update/independenttrials_03Nov2021_t012932/' ];
N_trials.mb_0p33( end + 1 )  = 10;

path_data.mb_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p0000001_oldsample_for_probe_update/independenttrials_02Nov2021_t120054/' ];
N_trials.mb_0p33( end + 1 )  = 10;

path_data.mb_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p00000001_oldsample_for_probe_update/independenttrials_01Nov2021_t234215/' ];
N_trials.mb_0p33( end + 1 )  = 10;

path_data.mb_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p000000001_oldsample_for_probe_update/independenttrials_01Nov2021_t113542/' ];
N_trials.mb_0p33( end + 1 )  = 10;

path_data.mb_0p33 = transpose( path_data.mb_0p33 );
N_trials.mb_0p33  = transpose( N_trials.mb_0p33 );

%==============
% 50% Minibatch
%==============

path_data.mb_0p50 = {};
N_trials.mb_0p50  = [];

path_data.mb_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha1p00_oldsample_for_probe_update/independenttrials_28Oct2021_t074706/' ];
N_trials.mb_0p50( end + 1 )  = 10;

path_data.mb_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p10_oldsample_for_probe_update/independenttrials_12Oct2021_t005707/' ];
N_trials.mb_0p50( end + 1 )  = 10;

path_data.mb_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p01_oldsample_for_probe_update/independenttrials_24Sep2021_t122324/' ];
N_trials.mb_0p50( end + 1 )  = 10;

path_data.mb_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p001_oldsample_for_probe_update/independenttrials_25Sep2021_t023411/' ];
N_trials.mb_0p50( end + 1 )  = 10;

path_data.mb_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p0001_oldsample_for_probe_update/independenttrials_25Sep2021_t164343/' ];
N_trials.mb_0p50( end + 1 )  = 10;

path_data.mb_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p00001_oldsample_for_probe_update/independenttrials_26Sep2021_t065757/' ];
N_trials.mb_0p50( end + 1 )  = 10;

path_data.mb_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p000001_oldsample_for_probe_update/independenttrials_02Nov2021_t102304/' ];
N_trials.mb_0p50( end + 1 )  = 10;

path_data.mb_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p0000001_oldsample_for_probe_update/independenttrials_03Nov2021_t031733/' ];
N_trials.mb_0p50( end + 1 )  = 10;

path_data.mb_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p00000001_oldsample_for_probe_update/independenttrials_03Nov2021_t201030/' ];
N_trials.mb_0p50( end + 1 )  = 10;

path_data.mb_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p000000001_oldsample_for_probe_update/independenttrials_04Nov2021_t125413/' ];
N_trials.mb_0p50( end + 1 )  = 10;

path_data.mb_0p50 = transpose( path_data.mb_0p50 );
N_trials.mb_0p50  = transpose( N_trials.mb_0p50 );

%==============
% Full Batch GD
%==============

path_data.full_gd = {};
N_trials.full_gd  = [];

path_data.full_gd{ end + 1 } = [ rootpath_data, 'full/cdi_rPIE_full_alpha1p00_randT_oldsample_for_probe_update/independenttrials_21Oct2021_t022747/' ];
N_trials.full_gd( end + 1 )  = 10;

path_data.full_gd{ end + 1 } = [ rootpath_data, 'full/cdi_rPIE_full_alpha0p10_randT_oldsample_for_probe_update/independenttrials_02Oct2021_t235705/' ];
N_trials.full_gd( end + 1 )  = 10;

path_data.full_gd{ end + 1 } = [ rootpath_data, 'full/cdi_rPIE_full_alpha0p01_randT_oldsample_for_probe_update/independenttrials_20Sep2021_t142149/' ];
N_trials.full_gd( end + 1 )  = 10;

path_data.full_gd{ end + 1 } = [ rootpath_data, 'full/cdi_rPIE_full_alpha0p001_randT_oldsample_for_probe_update/independenttrials_02Oct2021_t071559/' ];
N_trials.full_gd( end + 1 )  = 10;

path_data.full_gd{ end + 1 } = [ rootpath_data, 'full/cdi_rPIE_full_alpha0p0001_randT_oldsample_for_probe_update/independenttrials_10Nov2021_t120454/' ];
N_trials.full_gd( end + 1 )  = 10;

path_data.full_gd{ end + 1 } = [ rootpath_data, 'full/cdi_rPIE_full_alpha0p00001_randT_oldsample_for_probe_update/independenttrials_01Oct2021_t143516/' ];
N_trials.full_gd( end + 1 )  = 10;

path_data.full_gd{ end + 1 } = [ rootpath_data, 'full/cdi_rPIE_full_alpha0p000001_randT_oldsample_for_probe_update/independenttrials_11Nov2021_t080651/' ];
N_trials.full_gd( end + 1 )  = 10;

path_data.full_gd{ end + 1 } = [ rootpath_data, 'full/cdi_rPIE_full_alpha0p0000001_randT_oldsample_for_probe_update/independenttrials_30Sep2021_t215726/' ];
N_trials.full_gd( end + 1 )  = 10;

path_data.full_gd{ end + 1 } = [ rootpath_data, 'full/cdi_rPIE_full_alpha0p00000001_randT_oldsample_for_probe_update/independenttrials_12Nov2021_t023114/'];
N_trials.full_gd( end + 1 )  = 10;

path_data.full_gd{ end + 1 } = [ rootpath_data, 'full/cdi_rPIE_full_alpha0p000000001_randT_oldsample_for_probe_update/independenttrials_30Sep2021_t051649/' ];
N_trials.full_gd( end + 1 )  = 10;

path_data.full_gd = transpose( path_data.full_gd );
N_trials.full_gd  = transpose( N_trials.full_gd );

end

