%

%{

        cd /net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi

        rng( 'shuffle' )

        restoredefaultpath; 
        addpath( genpath( pwd ));
        clearvars -except expt sol




clear; close all; testing_mb_vs_rPIEalpha

%}

%=============================================================
% load previously processed metrics variable, do some plotting
%=============================================================

%{

Z = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/metrics_no_noise_MB_vs_rPIEalpha.mat' );

metrics       = Z.metrics;
probe_scaling = Z.probe_scaling;
plot_data_IDs = Z.path_data_fields;
N_data        = Z.N_data;
N_trials      = Z.N_trials;
path_data     = Z.path_data;

clear( 'Z' )

%========

N_alphatrials   = 10;
N_metricsepochs = 51;

metrics_plotting( length( plot_data_IDs ) ) = struct;

for jj = 1 : length( plot_data_IDs )                                    % loop over the 8 mini and full batch tests ( axis X )
    
%     metrics_plotting( jj ).data_id = ( plot_data_IDs{ jj } );

    for kk = 1 : length( path_data.( plot_data_IDs{ jj } ) )            % loop over the 8 MiniB + 2 FullB = 10 trials  
        
%         metrics_plotting( jj ).rPIEalpha{ kk } = metrics.( plot_data_IDs{ jj } ){ kk }.rPIE_alpha( 1 );       % 10 alpha test vals, { 1.0, 1e-1, ..., 1e-9 }
%     
%         metrics_plotting( jj ).log10_meas_all_avg{ kk } = metrics.( plot_data_IDs{ jj } ){ kk }.log10_meas_all_avg;
%         metrics_plotting( jj ).log10_meas_all_max{ kk } = transpose( metrics.( plot_data_IDs{ jj } ){ kk }.log10_meas_all_max );
%         metrics_plotting( jj ).log10_meas_all_min{ kk } = transpose( metrics.( plot_data_IDs{ jj } ){ kk }.log10_meas_all_min );
%         metrics_plotting( jj ).log10_meas_all_std{ kk } = std( metrics.( plot_data_IDs{ jj } ){ kk }.log10_meas_all, [], 2 );
%         
%         metrics_plotting( jj ).meas_all{ kk }           = metrics.( plot_data_IDs{ jj } ){ 2 }.meas_all;

        meas_all_avg( kk, jj, : ) = metrics.( plot_data_IDs{ jj } ){ kk }.meas_all_avg;
        meas_all_max( kk, jj, : ) = metrics.( plot_data_IDs{ jj } ){ kk }.meas_all_max;
        meas_all_min( kk, jj, : ) = metrics.( plot_data_IDs{ jj } ){ kk }.meas_all_min; 
        meas_all_std( kk, jj, : ) = std( metrics.( plot_data_IDs{ jj } ){ kk }.meas_all, [], 2 );
        
        log10_meas_all_avg( kk, jj, : ) = metrics.( plot_data_IDs{ jj } ){ kk }.log10_meas_all_avg;
        log10_meas_all_max( kk, jj, : ) = transpose( metrics.( plot_data_IDs{ jj } ){ kk }.log10_meas_all_max );
        log10_meas_all_min( kk, jj, : ) = transpose( metrics.( plot_data_IDs{ jj } ){ kk }.log10_meas_all_min );
        log10_meas_all_std( kk, jj, : ) = std( metrics.( plot_data_IDs{ jj } ){ kk }.log10_meas_all, [], 2 );
       
                
        for ii = 1 : 10
            
            timings_exwv_update( kk, jj, ii, : )   = metrics.( plot_data_IDs{ jj } ){ kk }.timing( ii ).exwv_update;
            timings_sample_update( kk, jj, ii, : ) = metrics.( plot_data_IDs{ jj } ){ kk }.timing( ii ).sample_update;
            timings_probe_update( kk, jj, ii, : )  = metrics.( plot_data_IDs{ jj } ){ kk }.timing( ii ).probe_update;
            timings_epoch( kk, jj, ii, : )         = metrics.( plot_data_IDs{ jj } ){ kk }.timing( ii ).epoch;
            
        end
        
        
        
        
    end

end

Y_alpha = [ 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9 ];
X_label = {'stoch\_bgd', 'mb\_0p01', 'mb\_0p05', 'mb\_0p10', 'mb\_0p20', 'mb\_0p33', 'mb\_0p50', 'full\_gd' }.';

clear( 'jj', 'ii', 'kk' )




% figure; 
% imagesc( log10( meas_all_avg( :, :, end ) ))
% xticklabels( X_label )
% yticklabels( Y_alpha )
% ylabel( 'rPIE alpha' )
% colormap turbo
% colorbar
% title('log10\_meas\_all\_avg')
% 
% figure; 
% imagesc( log10( meas_all_std( :, :, end ) ))
% xticklabels( X_label )
% yticklabels( Y_alpha )
% ylabel( 'rPIE alpha' )
% colormap turbo
% colorbar
% title('log10\_meas\_all\_std')
% 
% figure; 
% imagesc( log10( meas_all_max( :, :, end ) ))
% xticklabels( X_label )
% yticklabels( Y_alpha )
% ylabel( 'rPIE alpha' )
% colormap turbo
% colorbar
% title('log10\_meas\_all\_max')
% 
% figure; 
% imagesc( log10( meas_all_min( :, :, end ) ))
% xticklabels( X_label )
% yticklabels( Y_alpha )
% ylabel( 'rPIE alpha' )
% colormap turbo
% colorbar
% title('log10\_meas\_all\_min')






figure; 
imagesc( log10_meas_all_avg( :, :, end ) )
xticklabels( X_label )
yticklabels( Y_alpha )
ylabel( 'rPIE alpha' )
colormap turbo
colorbar
title('log10\_meas\_all\_avg')

figure; 
imagesc( log10_meas_all_std( :, :, end ) )
xticklabels( X_label )
yticklabels( Y_alpha )
ylabel( 'rPIE alpha' )
colormap turbo
colorbar
title('log10\_meas\_all\_std')

figure; 
imagesc( log10_meas_all_max( :, :, end ) )
xticklabels( X_label )
yticklabels( Y_alpha )
ylabel( 'rPIE alpha' )
colormap turbo
colorbar
title('log10\_meas\_all\_max')

figure; 
imagesc( log10_meas_all_min( :, :, end ) )
xticklabels( X_label )
yticklabels( Y_alpha )
ylabel( 'rPIE alpha' )
colormap turbo
colorbar
title('log10\_meas\_all\_min')


return
%}




%====================================================================================================================================================

% rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/';

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/noise_study/bg_rm/';
% rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/noise_study/no_bg_rm/';

%=============================================
% define the paths for the raw data to analyze
%=============================================

[ path_data, N_trials, N_data ] = define_paths( rootpath_data );

path_data_fields = fieldnames( path_data );

probe_scaling = ones( N_data, 1, 'single');
% probe_scaling = single( [ 1e-1, 1e0, 1e1, 1e2, 1e3 ] );
    

%====================================================================
% load and process the raw results, save to separate metrics variable
%====================================================================

for jj = 1 : length( path_data_fields )
 
    for kk = 1 : length( path_data.( path_data_fields{ jj } ) )
        
        metrics.( path_data_fields{ jj } ){ kk } = load_and_plot( path_data.( path_data_fields{ jj } ){ kk }, ...
                                                                  N_trials.( path_data_fields{ jj } )( kk ),  ...
                                                                  path_data_fields{ jj },                     ... 
                                                                  probe_scaling( kk )                           );                       
    
    end

end

%====================================================================================================================================================

function metrics = load_and_plot( path_data, N_trials, path_data_fields, probe_scaling )

    %========
    
    sim_ptycho2DTPA = cell( N_trials, 1 );
    
    for ii = 1 : N_trials

        sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
        
%         tmp0 = num2str( ii, 't%d' );
%         mkdir( tmp0 )
%         sim_ptycho2DTPA{ ii }.sol.print_img = true;
%         ptycho2DTPA_plotresults( sim_ptycho2DTPA{ ii }.sol, sim_ptycho2DTPA{ ii }.expt );
%         movefile( 'probe_5001.jpg', tmp0 )
%         movefile( 'sample_5001.jpg', tmp0 )

    end

    Nmetrics_epoch = length( sim_ptycho2DTPA{ 1 }.sol.metrics.meas_all );

    %========

    log10_meas_all_avg = 0;
    meas_all_avg       = 0;

    meas_all       = zeros( Nmetrics_epoch, N_trials, 'single' );
    log10_meas_all = zeros( Nmetrics_epoch, N_trials, 'single' );
    
    for ii = 1 : N_trials
        
        meas_all( :, ii ) = sim_ptycho2DTPA{ ii }.sol.metrics.meas_all / probe_scaling;     
%         meas_all( :, ii ) = sim_ptycho2DTPA{ ii }.sol.metrics.meas_all;     
        
        log10_meas_all( :, ii ) = log10( meas_all( :, ii ) ); 
        
        meas_all_avg       = meas_all_avg       + meas_all( :, ii );    
        log10_meas_all_avg = log10_meas_all_avg + log10_meas_all( :, ii );

    end

    log10_meas_all_avg = log10_meas_all_avg / N_trials;
    meas_all_avg       = meas_all_avg       / N_trials;
    
    log10_meas_all_min = transpose( min( log10_meas_all, [], 2 ));
    log10_meas_all_max = transpose( max( log10_meas_all, [], 2 ));  
    
    meas_all_min = transpose( min( meas_all, [], 2 ));
    meas_all_max = transpose( max( meas_all, [], 2 ));  
    
%     name_data = num2str( [ sim_ptycho2DTPA{ 1 }.sol.rPIE_alpha, sim_ptycho2DTPA{ 1 }.sol.spos.rand_spos_subset_pct ], 'rPIE_alpha = %0.8f, MBpct = %0.4f');

    %=================================================
    % create a struct from metrics data for future use
    %=================================================
    
    for ii = 1 : N_trials
        
        metrics.meas_all( :, ii )     = sim_ptycho2DTPA{ ii }.sol.metrics.meas_all;
        metrics.timing( ii )          = sim_ptycho2DTPA{ ii }.sol.timings;
        metrics.it( ii )              = sim_ptycho2DTPA{ ii }.sol.it;
        metrics.rPIE_alpha( ii )      = sim_ptycho2DTPA{ ii }.sol.rPIE_alpha;
        metrics.probe_intensity( ii ) = sim_ptycho2DTPA{ ii }.sol.probe.scpm.fro2TOT;
        
        if isfield( sim_ptycho2DTPA{ ii }.sol.spos, 'rand_spos_subset_pct' )
            
            metrics.rand_spos_subset_pct( ii ) = sim_ptycho2DTPA{ ii }.sol.spos.rand_spos_subset_pct;
            
        end
        
    end
    
    %========
    
    if isfield( sim_ptycho2DTPA{ ii }.sol.spos, 'rand_spos_subset_pct' )
        
        name_data = [ path_data_fields, num2str( [ metrics.rPIE_alpha( 1 ),             ...
                                                   metrics.rand_spos_subset_pct( 1 ) ,  ...
                                                   round( metrics.probe_intensity( 1 )) ], ', rPIE_alpha = %0.9f, MBpct = %0.4f, probe photons = %0.3e') ];
                                               
        save_data = [ path_data_fields, num2str( [ metrics.rPIE_alpha( 1 ),             ...
                                                   metrics.rand_spos_subset_pct( 1 ) ,  ...
                                                   round( metrics.probe_intensity( 1 )) ], '_rPIEalpha_%0.9f_MBpct_%0.4f_photons_%0.3e'), '.jpg' ];
        
    else
        
        name_data = [ path_data_fields, num2str( [ sim_ptycho2DTPA{ ii }.sol.rPIE_alpha, ...
                                                   round( metrics.probe_intensity( 1 )) ], ', rPIE_alpha = %0.9f, probe photons = %0.3e') ];
                                               
        save_data = [ path_data_fields, num2str( [ sim_ptycho2DTPA{ ii }.sol.rPIE_alpha, ...
                                                   round( metrics.probe_intensity( 1 )) ], '_rPIEalpha_%0.9f_photons_%0.3e'), '.jpg' ];
        
    end
    
    %========
    
    metrics.log10_meas_all_max  = log10_meas_all_max;
    metrics.log10_meas_all_min  = log10_meas_all_min;
    
    metrics.meas_all_max        = meas_all_max;
    metrics.meas_all_min        = meas_all_min;
    
    metrics.meas_all            = meas_all;
    metrics.log10_meas_all      = log10_meas_all;
    
    metrics.log10_meas_all_avg  = log10_meas_all_avg;
    metrics.meas_all_avg        = meas_all_avg;
    
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
%     plot( metrics.it( 1 ).mtot( 1 : skip : end ), log10_meas_all_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
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
     
    %==========================================================================
    % normalize metrics vs epoch trajectory plots so that starting error is 1.0
    %==========================================================================

    normalize_metric_start = logical( 0 ); %#ok<LOGL>

    if normalize_metric_start == true
        
        tmp0 = log10( metrics.meas_all / probe_scaling );
        tmp1 = log10_meas_all_avg      / probe_scaling;
 
    else

        tmp0 = log10( metrics.meas_all );
        tmp1 = log10_meas_all_avg;
    end

    tmp2 = metrics.it( 1 ).mtot;
    
    %=========================
    % create the plot, save it
    %=========================
    
    skip = 1;
    
    y_lim = [ -2, 9 ];
%     y_lim = [ 4, 7 ];

    h1 = figure();  
    set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )

    hold on

    for ii = 1 : N_trials

        plot( tmp2( 1 : skip : end ), tmp0( 1 : skip : end, ii ), '-', 'linewidth', 2 )

    end

    plot( tmp2( 1 : skip : end ), tmp1( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )

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

function [ path_data, N_trials, N_data ] = define_paths( rootpath_data )

%===============
% initialization
%===============

N_data = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               Noisy vs Non-noisy measurements 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%===============
% NO NOISE MB 1%
%===============

%===============
% NO NOISE MB 5%
%===============

%================
% NO NOISE MB 10%
%================
% 
% path_data.mb_0p01 = {};
% N_trials.mb_0p01  = [];
% 
% path_data.mb_0p01{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p01_oldsample_for_probe_update/independenttrials_20Sep2021_t143405/' ];
% N_trials.mb_0p01( end + 1 )  = 10;
% 
% path_data.mb_0p01{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p001_oldsample_for_probe_update/independenttrials_20Sep2021_t144026/' ];
% N_trials.mb_0p01( end + 1 )  = 10;
% 
% path_data.mb_0p01 = transpose( path_data.mb_0p01 );
% N_trials.mb_0p01  = transpose( N_trials.mb_0p01 );

% N_data = N_data + length( path_data.mb_0p01 );

%================
% NO NOISE MB 20%
%================

%================
% NO NOISE MB 33%
%================

%================
% NO NOISE MB 50%
%================

%=================
% NO NOISE FULL GD
%=================

%%%%%%%%%%%%%%%%%%% WITH NOISE, NO BG RM %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% WITH NOISE, NO BG RM %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% WITH NOISE, NO BG RM %%%%%%%%%%%%%%%%%%%

% %=====================
% % NOISY Block Stoch GD
% %=====================
% 
% path_data.blockstochGD_noise = {};
% N_trials.blockstochGD_noise  = [];
% 
% path_data.blockstochGD_noise{ end + 1 } = [ rootpath_data, 'cdi_rPIE_stochGD_alpha0p10_noise/independenttrials_19Nov2021_t072323/' ];
% N_trials.blockstochGD_noise( end + 1 )  = 10;
% 
% path_data.blockstochGD_noise{ end + 1 } = [ rootpath_data, 'cdi_rPIE_stochGD_alpha0p01_noise/independenttrials_17Nov2021_t152851/' ];
% N_trials.blockstochGD_noise( end + 1 )  = 10;
% 
% path_data.blockstochGD_noise{ end + 1 } = [ rootpath_data, 'cdi_rPIE_stochGD_alpha0p001_noise/independenttrials_18Nov2021_t112543/' ];
% N_trials.blockstochGD_noise( end + 1 )  = 10;
% 
% path_data.blockstochGD_noise = transpose( path_data.blockstochGD_noise );
% N_trials.blockstochGD_noise  = transpose( N_trials.blockstochGD_noise );
% 
% N_data = N_data + length( path_data.blockstochGD_noise );
% 
% %============
% % NOISY MB 1%
% %============
% 
% %============
% % NOISY MB 5%
% %============
% 
% %=============
% % NOISY MB 10%
% %=============
% 
% path_data.mb_0p10_noise = {};
% N_trials.mb_0p10_noise  = [];
% 
% path_data.mb_0p10_noise{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p10_alpha1p0_noise/independenttrials_16Nov2021_t173234/' ];
% N_trials.mb_0p10_noise( end + 1 )  = 10;
% 
% path_data.mb_0p10_noise{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p10_alpha0p1_noise/independenttrials_16Nov2021_t091233/' ];
% N_trials.mb_0p10_noise( end + 1 )  = 10;
% 
% path_data.mb_0p10_noise{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p10_alpha0p01_noise/independenttrials_15Nov2021_t144212/' ];
% N_trials.mb_0p10_noise( end + 1 )  = 10;
% 
% path_data.mb_0p10_noise{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p10_alpha0p001_noise/independenttrials_15Nov2021_t230419/' ];
% N_trials.mb_0p10_noise( end + 1 )  = 10;
% 
% path_data.mb_0p10_noise = transpose( path_data.mb_0p10_noise );
% N_trials.mb_0p10_noise  = transpose( N_trials.mb_0p10_noise );
% 
% N_data = N_data + length( path_data.mb_0p10_noise );
% 
% %=============
% % NOISY MB 20%
% %=============
% 
% path_data.mb_0p20_noise = {};
% N_trials.mb_0p20_noise  = [];
% 
% path_data.mb_0p20_noise{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p20_alpha1p0_noise/independenttrials_16Nov2021_t190023/' ];
% N_trials.mb_0p20_noise( end + 1 )  = 10;
% 
% path_data.mb_0p20_noise{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p20_alpha0p1_noise/independenttrials_16Nov2021_t091622/' ];
% N_trials.mb_0p20_noise( end + 1 )  = 10;
% 
% path_data.mb_0p20_noise = transpose( path_data.mb_0p20_noise );
% N_trials.mb_0p20_noise  = transpose( N_trials.mb_0p20_noise );
% 
% N_data = N_data + length( path_data.mb_0p20_noise );
% 
% %=============
% % NOISY MB 33%
% %=============
% 
% path_data.mb_0p33_noise = {};
% N_trials.mb_0p33_noise  = [];
% 
% path_data.mb_0p33_noise{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p33_alpha1p0_noise/independenttrials_16Nov2021_t193743/' ];
% N_trials.mb_0p33_noise( end + 1 )  = 10;
% 
% path_data.mb_0p33_noise{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p33_alpha0p1_noise/independenttrials_16Nov2021_t091946/' ];
% N_trials.mb_0p33_noise( end + 1 )  = 10;
% 
% path_data.mb_0p33_noise = transpose( path_data.mb_0p33_noise );
% N_trials.mb_0p33_noise  = transpose( N_trials.mb_0p33_noise );
% 
% N_data = N_data + length( path_data.mb_0p33_noise );
% 
% %=============
% % NOISY MB 50%
% %=============
% 
% path_data.mb_0p50_noise = {};
% N_trials.mb_0p50_noise  = [];
% 
% path_data.mb_0p50_noise{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p50_alpha0p01_noise/independenttrials_17Nov2021_t152800/' ];
% N_trials.mb_0p50_noise( end + 1 )  = 10;
% 
% path_data.mb_0p50_noise{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p50_alpha0p001_noise/independenttrials_18Nov2021_t153057/' ];
% N_trials.mb_0p50_noise( end + 1 )  = 10;
% 
% path_data.mb_0p50_noise = transpose( path_data.mb_0p50_noise );
% N_trials.mb_0p50_noise  = transpose( N_trials.mb_0p50_noise );
% 
% N_data = N_data + length( path_data.mb_0p50_noise );
% 
% %==============
% % NOISY FULL GD
% %==============
% 
% path_data.fullGD_noise = {};
% N_trials.fullGD_noise  = [];
% 
% path_data.fullGD_noise{ end + 1 } = [ rootpath_data, 'cdi_rPIE_full_alpha1p00_randT_noise/independenttrials_19Nov2021_t193738/' ];
% N_trials.fullGD_noise( end + 1 )  = 10;
% 
% path_data.fullGD_noise{ end + 1 } = [ rootpath_data, 'cdi_rPIE_full_alpha0p10_randT_noise/independenttrials_19Nov2021_t021441/' ];
% N_trials.fullGD_noise( end + 1 )  = 10;
% 
% path_data.fullGD_noise{ end + 1 } = [ rootpath_data, 'cdi_rPIE_full_alpha0p01_randT_noise/independenttrials_17Nov2021_t152829/' ];
% N_trials.fullGD_noise( end + 1 )  = 10;
% 
% path_data.fullGD_noise{ end + 1 } = [ rootpath_data, 'cdi_rPIE_full_alpha0p001_randT_noise/independenttrials_18Nov2021_t085030/' ];
% N_trials.fullGD_noise( end + 1 )  = 10;
% 
% path_data.fullGD_noise = transpose( path_data.fullGD_noise );
% N_trials.fullGD_noise  = transpose( N_trials.fullGD_noise );
% 
% N_data = N_data + length( path_data.fullGD_noise );
% 
% 
% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WITH NOISE, WITH BG RM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WITH NOISE, WITH BG RM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WITH NOISE, WITH BG RM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%============
% NOISY MB 1%
%============

% path_data.mb_0p01_noise_rmbg = {};
% N_trials.mb_0p01_noise_rmbg  = [];
% 
% % path_data.mb_0p01_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p01_alpha1p00_noise_rmbg/independenttrials_24Nov2021_t170431/' ];
% % N_trials.mb_0p01_noise_rmbg( end + 1 )  = 10;
% % 
% % path_data.mb_0p01_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p01_alpha0p10_noise_rmbg/independenttrials_26Nov2021_t073554/' ];
% % N_trials.mb_0p01_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p01_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p01_alpha0p01_noise_rmbg/independenttrials_27Nov2021_t220437/' ];
% N_trials.mb_0p01_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p01_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p01_alpha0p001_noise_rmbg/independenttrials_29Nov2021_t122721/' ];
% N_trials.mb_0p01_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p01_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p01_alpha0p0001_noise_rmbg/independenttrials_01Dec2021_t010426/' ];
% N_trials.mb_0p01_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p01_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p01_alpha0p00001_noise_rmbg/independenttrials_01Dec2021_t231204/' ];
% N_trials.mb_0p01_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p01_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p01_alpha0p000001_noise_rmbg/independenttrials_02Dec2021_t213239/' ];
% N_trials.mb_0p01_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p01_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p01_alpha0p0000001_noise_rmbg/independenttrials_03Dec2021_t204057/' ];
% N_trials.mb_0p01_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p01_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p01_alpha0p00000001_noise_rmbg/independenttrials_04Dec2021_t184438/' ];
% N_trials.mb_0p01_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p01_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p01_alpha0p000000001_noise_rmbg/independenttrials_05Dec2021_t164611/' ];
% N_trials.mb_0p01_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p01_noise_rmbg = transpose( path_data.mb_0p01_noise_rmbg );
% N_trials.mb_0p01_noise_rmbg  = transpose( N_trials.mb_0p01_noise_rmbg );
% 
% N_data = N_data + length( path_data.mb_0p01_noise_rmbg );

%============================================
% NOISY MB 1% W/ GROUND TRUTH SOLUTION STARTS
%============================================

% path_data.mb_0p01_noise_rmbg_soln = {};
% N_trials.mb_0p01_noise_rmbg_soln  = [];
% 
% path_data.mb_0p01_noise_rmbg_soln{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p01_alpha0p000000001_noise_rmbg_SOLUTION/independenttrials_25Nov2021_t151935/' ];
% N_trials.mb_0p01_noise_rmbg_soln( end + 1 )  = 10;
% 
% path_data.mb_0p01_noise_rmbg_soln = transpose( path_data.mb_0p01_noise_rmbg_soln );
% N_trials.mb_0p01_noise_rmbg_soln  = transpose( N_trials.mb_0p01_noise_rmbg_soln );
% 
% N_data = N_data + length( path_data.mb_0p01_noise_rmbg_soln );

%============
% NOISY MB 5%
%============

% path_data.mb_0p05_noise_rmbg = {};
% N_trials.mb_0p05_noise_rmbg  = [];
% 
% path_data.mb_0p05_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p05_alpha1p00_noise_rmbg/independenttrials_01Dec2021_t192641/' ];
% N_trials.mb_0p05_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p05_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p05_alpha0p10_noise_rmbg/independenttrials_02Dec2021_t065029/' ];
% N_trials.mb_0p05_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p05_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p05_alpha0p01_noise_rmbg/independenttrials_06Dec2021_t020602/' ];
% N_trials.mb_0p05_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p05_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p05_alpha0p001_noise_rmbg/independenttrials_05Dec2021_t145000/' ];
% N_trials.mb_0p05_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p05_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p05_alpha0p0001_noise_rmbg/independenttrials_05Dec2021_t033413/' ];
% N_trials.mb_0p05_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p05_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p05_alpha0p00001_noise_rmbg/independenttrials_04Dec2021_t161728/' ];
% N_trials.mb_0p05_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p05_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p05_alpha0p000001_noise_rmbg/independenttrials_04Dec2021_t045955/' ];
% N_trials.mb_0p05_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p05_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p05_alpha0p0000001_noise_rmbg/independenttrials_03Dec2021_t174213/' ];
% N_trials.mb_0p05_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p05_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p05_alpha0p00000001_noise_rmbg/independenttrials_03Dec2021_t055501/' ];
% N_trials.mb_0p05_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p05_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p05_alpha0p000000001_noise_rmbg/independenttrials_02Dec2021_t181433/' ];
% N_trials.mb_0p05_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p05_noise_rmbg = transpose( path_data.mb_0p05_noise_rmbg );
% N_trials.mb_0p05_noise_rmbg  = transpose( N_trials.mb_0p05_noise_rmbg );
% 
% N_data = N_data + length( path_data.mb_0p05_noise_rmbg );
% 
% return

%============================================
% NOISY MB 5% W/ GROUND TRUTH SOLUTION STARTS
%============================================

% path_data.mb_0p05_noise_rmbg_soln = {};
% N_trials.mb_0p05_noise_rmbg_soln  = [];
% 
% path_data.mb_0p05_noise_rmbg_soln{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p05_alpha0p000000001_noise_rmbg_SOLUTION/independenttrials_27Nov2021_t011837/' ];
% N_trials.mb_0p05_noise_rmbg_soln( end + 1 )  = 10;
% 
% path_data.mb_0p05_noise_rmbg_soln = transpose( path_data.mb_0p05_noise_rmbg_soln );
% N_trials.mb_0p05_noise_rmbg_soln  = transpose( N_trials.mb_0p05_noise_rmbg_soln );
% 
% N_data = N_data + length( path_data.mb_0p05_noise_rmbg_soln );

%=============
% NOISY MB 10%
%=============

% path_data.mb_0p10_noise_rmbg = {};
% N_trials.mb_0p10_noise_rmbg  = [];
% 
% path_data.mb_0p10_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p10_alpha1p00_noise_rmbg/independenttrials_22Nov2021_t004230/' ];
% N_trials.mb_0p10_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p10_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p10_alpha0p10_noise_rmbg/independenttrials_21Nov2021_t161037/' ];
% N_trials.mb_0p10_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p10_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p10_alpha0p01_noise_rmbg/independenttrials_20Nov2021_t154800/' ];
% N_trials.mb_0p10_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p10_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p10_alpha0p001_noise_rmbg/independenttrials_21Nov2021_t040200/' ];
% N_trials.mb_0p10_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p10_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p10_alpha0p0001_noise_rmbg/independenttrials_01Dec2021_t092544/' ];
% N_trials.mb_0p10_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p10_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p10_alpha0p00001_noise_rmbg/independenttrials_30Nov2021_t233736/' ];
% N_trials.mb_0p10_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p10_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p10_alpha0p000001_noise_rmbg/independenttrials_30Nov2021_t085022/' ];
% N_trials.mb_0p10_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p10_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p10_alpha0p0000001_noise_rmbg/independenttrials_29Nov2021_t162031/' ];
% N_trials.mb_0p10_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p10_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p10_alpha0p00000001_noise_rmbg/independenttrials_28Nov2021_t235054/' ];
% N_trials.mb_0p10_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p10_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p10_alpha0p000000001_noise_rmbg/independenttrials_28Nov2021_t072021/' ];
% N_trials.mb_0p10_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p10_noise_rmbg = transpose( path_data.mb_0p10_noise_rmbg );
% N_trials.mb_0p10_noise_rmbg  = transpose( N_trials.mb_0p10_noise_rmbg );
% 
% N_data = N_data + length( path_data.mb_0p10_noise_rmbg );

%=============================================
% NOISY MB 10% W/ GROUND TRUTH SOLUTION STARTS
%=============================================

% path_data.mb_0p10_noise_rmbg_soln = {};
% N_trials.mb_0p10_noise_rmbg_soln  = [];
% 
% path_data.mb_0p10_noise_rmbg_soln{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p10_alpha0p000000001_noise_rmbg_SOLUTION/independenttrials_22Nov2021_t223444/' ];
% N_trials.mb_0p10_noise_rmbg_soln( end + 1 )  = 10;
% 
% path_data.mb_0p10_noise_rmbg_soln = transpose( path_data.mb_0p10_noise_rmbg_soln );
% N_trials.mb_0p10_noise_rmbg_soln  = transpose( N_trials.mb_0p10_noise_rmbg_soln );
% 
% N_data = N_data + length( path_data.mb_0p10_noise_rmbg_soln );

%=============
% NOISY MB 20%
%=============

% path_data.mb_0p20_noise_rmbg = {};
% N_trials.mb_0p20_noise_rmbg  = [];
% 
% path_data.mb_0p20_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p20_alpha1p00_noise_rmbg/independenttrials_22Nov2021_t054607/' ];
% N_trials.mb_0p20_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p20_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p20_alpha0p10_noise_rmbg/independenttrials_21Nov2021_t195559/' ];
% N_trials.mb_0p20_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p20_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p20_alpha0p01_noise_rmbg/independenttrials_20Nov2021_t160330/' ];
% N_trials.mb_0p20_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p20_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p20_alpha0p001_noise_rmbg/independenttrials_21Nov2021_t075000/' ];
% N_trials.mb_0p20_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p20_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p20_alpha0p0001_noise_rmbg/independenttrials_27Nov2021_t124329/' ];
% N_trials.mb_0p20_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p20_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p20_alpha0p00001_noise_rmbg/independenttrials_26Nov2021_t180859/' ];
% N_trials.mb_0p20_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p20_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p20_alpha0p000001_noise_rmbg/independenttrials_25Nov2021_t233530/' ];
% N_trials.mb_0p20_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p20_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p20_alpha0p0000001_noise_rmbg/independenttrials_25Nov2021_t050055/' ];
% N_trials.mb_0p20_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p20_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p20_alpha0p00000001_noise_rmbg/independenttrials_24Nov2021_t102832/' ];
% N_trials.mb_0p20_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p20_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p20_alpha0p000000001_noise_rmbg/independenttrials_23Nov2021_t160745/' ];
% N_trials.mb_0p20_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p20_noise_rmbg = transpose( path_data.mb_0p20_noise_rmbg );
% N_trials.mb_0p20_noise_rmbg  = transpose( N_trials.mb_0p20_noise_rmbg );
% 
% N_data = N_data + length( path_data.mb_0p20_noise_rmbg );

%=============================================
% NOISY MB 20% W/ GROUND TRUTH SOLUTION STARTS
%=============================================

% path_data.mb_0p20_noise_rmbg_soln = {};
% N_trials.mb_0p20_noise_rmbg_soln  = [];
% 
% path_data.mb_0p20_noise_rmbg_soln{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p20_alpha0p000000001_noise_rmbg_SOLUTION/independenttrials_22Nov2021_t130028/' ];
% N_trials.mb_0p20_noise_rmbg_soln( end + 1 )  = 10;
% 
% path_data.mb_0p20_noise_rmbg_soln = transpose( path_data.mb_0p20_noise_rmbg_soln );
% N_trials.mb_0p20_noise_rmbg_soln  = transpose( N_trials.mb_0p20_noise_rmbg_soln );
% 
% N_data = N_data + length( path_data.mb_0p20_noise_rmbg_soln );

%=============
% NOISY MB 33%
%=============

% path_data.mb_0p33_noise_rmbg = {};
% N_trials.mb_0p33_noise_rmbg  = [];
% 
% path_data.mb_0p33_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p33_alpha1p00_noise_rmbg/independenttrials_22Nov2021_t120012/' ];
% N_trials.mb_0p33_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p33_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p33_alpha0p10_noise_rmbg/independenttrials_22Nov2021_t013955/' ];
% N_trials.mb_0p33_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p33_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p33_alpha0p01_noise_rmbg/independenttrials_20Nov2021_t160810/' ];
% N_trials.mb_0p33_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p33_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p33_alpha0p001_noise_rmbg/independenttrials_21Nov2021_t142752/' ];
% N_trials.mb_0p33_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p33_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p33_alpha0p0001_noise_rmbg/independenttrials_28Nov2021_t220333/' ];
% N_trials.mb_0p33_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p33_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p33_alpha0p00001_noise_rmbg/independenttrials_29Nov2021_t173018/' ];
% N_trials.mb_0p33_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p33_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p33_alpha0p000001_noise_rmbg/independenttrials_30Nov2021_t125707/' ];
% N_trials.mb_0p33_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p33_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p33_alpha0p0000001_noise_rmbg/independenttrials_01Dec2021_t040331/' ];
% N_trials.mb_0p33_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p33_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p33_alpha0p00000001_noise_rmbg/independenttrials_01Dec2021_t161052/' ];
% N_trials.mb_0p33_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p33_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p33_alpha0p000000001_noise_rmbg/independenttrials_02Dec2021_t042237/' ];
% N_trials.mb_0p33_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p33_noise_rmbg = transpose( path_data.mb_0p33_noise_rmbg );
% N_trials.mb_0p33_noise_rmbg  = transpose( N_trials.mb_0p33_noise_rmbg );
% 
% N_data = N_data + length( path_data.mb_0p33_noise_rmbg );

%=============================================
% NOISY MB 33% W/ GROUND TRUTH SOLUTION STARTS
%=============================================

% path_data.mb_0p33_noise_rmbg_soln = {};
% N_trials.mb_0p33_noise_rmbg_soln  = [];
% 
% path_data.mb_0p33_noise_rmbg_soln{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p33_alpha0p000000001_noise_rmbg_SOLUTION/independenttrials_02Dec2021_t163415/' ];
% N_trials.mb_0p33_noise_rmbg_soln( end + 1 )  = 10;
% 
% path_data.mb_0p33_noise_rmbg_soln = transpose( path_data.mb_0p33_noise_rmbg_soln );
% N_trials.mb_0p33_noise_rmbg_soln  = transpose( N_trials.mb_0p33_noise_rmbg_soln );
% 
% N_data = N_data + length( path_data.mb_0p33_noise_rmbg_soln );

%=============
% NOISY MB 50%
%=============

% path_data.mb_0p50_noise_rmbg = {};
% N_trials.mb_0p50_noise_rmbg  = [];
% 
% path_data.mb_0p50_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p50_alpha1p00_noise_rmbg/independenttrials_30Nov2021_t211459/' ];
% N_trials.mb_0p50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p50_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p50_alpha0p10_noise_rmbg/independenttrials_01Dec2021_t124737/' ];
% N_trials.mb_0p50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p50_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p50_alpha0p01_noise_rmbg/independenttrials_02Dec2021_t041855/' ];
% N_trials.mb_0p50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p50_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p50_alpha0p001_noise_rmbg/independenttrials_02Dec2021_t195009/' ];
% N_trials.mb_0p50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p50_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p50_alpha0p0001_noise_rmbg/independenttrials_03Dec2021_t114103/' ];
% N_trials.mb_0p50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p50_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p50_alpha0p00001_noise_rmbg/independenttrials_04Dec2021_t031455/' ];
% N_trials.mb_0p50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p50_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p50_alpha0p000001_noise_rmbg/independenttrials_04Dec2021_t184436/' ];
% N_trials.mb_0p50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p50_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p50_alpha0p0000001_noise_rmbg/independenttrials_05Dec2021_t101235/' ];
% N_trials.mb_0p50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p50_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p50_alpha0p00000001_noise_rmbg/independenttrials_06Dec2021_t014006/' ];
% N_trials.mb_0p50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p50_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p50_alpha0p000000001_noise_rmbg/independenttrials_06Dec2021_t170610/' ];
% N_trials.mb_0p50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb_0p50_noise_rmbg = transpose( path_data.mb_0p50_noise_rmbg );
% N_trials.mb_0p50_noise_rmbg  = transpose( N_trials.mb_0p50_noise_rmbg );
% 
% N_data = N_data + length( path_data.mb_0p50_noise_rmbg );

%=============================================
% NOISY MB 50% W/ GROUND TRUTH SOLUTION STARTS
%=============================================

% path_data.mb_0p50_noise_rmbg_soln = {};
% N_trials.mb_0p50_noise_rmbg_soln  = [];
% 
% path_data.mb_0p50_noise_rmbg_soln{ end + 1 } = [ rootpath_data, 'cdi_rPIE_mb0p50_alpha0p000000001_noise_rmbg_SOLUTION/independenttrials_03Dec2021_t003728/' ];
% N_trials.mb_0p50_noise_rmbg_soln( end + 1 )  = 10;
% 
% path_data.mb_0p50_noise_rmbg_soln = transpose( path_data.mb_0p50_noise_rmbg_soln );
% N_trials.mb_0p50_noise_rmbg_soln  = transpose( N_trials.mb_0p50_noise_rmbg_soln );
% 
% N_data = N_data + length( path_data.mb_0p50_noise_rmbg_soln );

%===================
% NOISY stochblockGD
%===================

% path_data.blockstochGD_noise_rmbg = {};
% N_trials.blockstochGD_noise_rmbg  = [];
% 
% path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_stochGD_alpha1p00_noise_rmbg/independenttrials_24Nov2021_t180523/' ];
% N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_stochGD_alpha0p10_noise_rmbg/independenttrials_23Nov2021_t213427/' ];
% N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_stochGD_alpha0p01_noise_rmbg/independenttrials_22Nov2021_t120907/' ];
% N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_stochGD_alpha0p001_noise_rmbg/independenttrials_23Nov2021_t041408/' ];
% N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_stochGD_alpha0p0001_noise_rmbg/independenttrials_03Dec2021_t135839/' ];
% N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_stochGD_alpha0p00001_noise_rmbg/independenttrials_04Dec2021_t070508/' ];
% N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_stochGD_alpha0p000001_noise_rmbg/independenttrials_05Dec2021_t001221/' ];
% N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_stochGD_alpha0p0000001_noise_rmbg/independenttrials_05Dec2021_t171931/' ];
% N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_stochGD_alpha0p00000001_noise_rmbg/independenttrials_06Dec2021_t102516/' ];
% N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_stochGD_alpha0p000000001_noise_rmbg/independenttrials_07Dec2021_t031749/' ];
% N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.blockstochGD_noise_rmbg = transpose( path_data.blockstochGD_noise_rmbg );
% N_trials.blockstochGD_noise_rmbg  = transpose( N_trials.blockstochGD_noise_rmbg );
% 
% N_data = N_data + length( path_data.blockstochGD_noise_rmbg );

%===================================================
% NOISY stochblockGD W/ GROUND TRUTH SOLUTION STARTS
%===================================================

% path_data.blockstochGD_noise_rmbg_soln = {};
% N_trials.blockstochGD_noise_rmbg_soln  = [];
% 
% path_data.blockstochGD_noise_rmbg_soln{ end + 1 } = [ rootpath_data, 'cdi_rPIE_stochGD_alpha0p000000001_noise_rmbg_SOLUTION/independenttrials_07Dec2021_t192851/' ];
% N_trials.blockstochGD_noise_rmbg_soln( end + 1 )  = 10;
% 
% path_data.blockstochGD_noise_rmbg_soln = transpose( path_data.blockstochGD_noise_rmbg_soln );
% N_trials.blockstochGD_noise_rmbg_soln  = transpose( N_trials.blockstochGD_noise_rmbg_soln );
% 
% N_data = N_data + length( path_data.blockstochGD_noise_rmbg_soln );

%==================
% NOISY fullbatchGD
%==================

% path_data.fullGD_noise_rmbg = {};
% N_trials.fullGD_noise_rmbg  = [];
% 
% path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_full_alpha1p00_randT_noise_rmbg/independenttrials_23Nov2021_t160603/' ];
% N_trials.fullGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_full_alpha0p10_randT_noise_rmbg/independenttrials_24Nov2021_t091610/' ];
% N_trials.fullGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_full_alpha0p01_randT_noise_rmbg/independenttrials_25Nov2021_t023340/' ];
% N_trials.fullGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_full_alpha0p001_randT_noise_rmbg/independenttrials_25Nov2021_t195309/' ];
% N_trials.fullGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_full_alpha0p0001_randT_noise_rmbg/independenttrials_26Nov2021_t131152/' ];
% N_trials.fullGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_full_alpha0p00001_randT_noise_rmbg/independenttrials_27Nov2021_t063024/' ];
% N_trials.fullGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_full_alpha0p000001_randT_noise_rmbg/independenttrials_27Nov2021_t235023/' ];
% N_trials.fullGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_full_alpha0p0000001_randT_noise_rmbg/independenttrials_28Nov2021_t171116/' ];
% N_trials.fullGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_full_alpha0p00000001_randT_noise_rmbg/independenttrials_29Nov2021_t103218/' ];
% N_trials.fullGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'cdi_rPIE_full_alpha0p000000001_randT_noise_rmbg/independenttrials_30Nov2021_t035342/' ];
% N_trials.fullGD_noise_rmbg( end + 1 )  = 10;
% 
% path_data.fullGD_noise_rmbg = transpose( path_data.fullGD_noise_rmbg );
% N_trials.fullGD_noise_rmbg  = transpose( N_trials.fullGD_noise_rmbg );
% 
% N_data = N_data + length( path_data.fullGD_noise_rmbg );

%==================================================
% NOISY fullbatchGD W/ GROUND TRUTH SOLUTION STARTS
%==================================================

% path_data.fullGD_noise_rmbg_soln = {};
% N_trials.fullGD_noise_rmbg_soln  = [];
% 
% path_data.fullGD_noise_rmbg_soln{ end + 1 } = [ rootpath_data, 'cdi_rPIE_full_alpha0p000000001_randT_noise_rmbg_SOLUTION/independenttrials_07Dec2021_t082247/' ];
% N_trials.fullGD_noise_rmbg_soln( end + 1 )  = 10;
% 
% path_data.fullGD_noise_rmbg_soln = transpose( path_data.fullGD_noise_rmbg_soln );
% N_trials.fullGD_noise_rmbg_soln  = transpose( N_trials.fullGD_noise_rmbg_soln );
% 
% N_data = N_data + length( path_data.fullGD_noise_rmbg_soln );

return


%==================================================================================================
%                                       Probe Scaling Study
%==================================================================================================

% path_data.mb_0p10_nonoise = {};
% N_trials.mb_0p10_nonoise  = [];
% 
% path_data.mb_0p10_nonoise{ end + 1 } = [ rootpath_data, '/scaling_study/nonoise/cdi_rPIE_mb0p10_alpha0p001_1e-1_probe_scaling/independenttrials_15Nov2021_t144109/' ];
% N_trials.mb_0p10_nonoise( end + 1 )  = 10;
% 
% path_data.mb_0p10_nonoise{ end + 1 } = [ rootpath_data, '/mb0p10/cdi_rPIE_mb0p10_alpha0p001_oldsample_for_probe_update/independenttrials_20Sep2021_t144026/' ];
% N_trials.mb_0p10_nonoise( end + 1 )  = 10;
% 
% path_data.mb_0p10_nonoise{ end + 1 } = [ rootpath_data, '/scaling_study/nonoise/cdi_rPIE_mb0p10_alpha0p001_1e+1_probe_scaling/independenttrials_15Nov2021_t231206/' ];
% N_trials.mb_0p10_nonoise( end + 1 )  = 10;
% 
% path_data.mb_0p10_nonoise{ end + 1 } = [ rootpath_data, '/scaling_study/nonoise/cdi_rPIE_mb0p10_alpha0p001_1e+2_probe_scaling/independenttrials_15Nov2021_t144304/' ];
% N_trials.mb_0p10_nonoise( end + 1 )  = 10;
% 
% path_data.mb_0p10_nonoise{ end + 1 } = [ rootpath_data, '/scaling_study/nonoise/cdi_rPIE_mb0p10_alpha0p001_1e+3_probe_scaling/independenttrials_15Nov2021_t230933/' ];
% N_trials.mb_0p10_nonoise( end + 1 )  = 10;
% 
% path_data.mb_0p10_nonoise = transpose( path_data.mb_0p10_nonoise );
% N_trials.mb_0p10_nonoise  = transpose( N_trials.mb_0p10_nonoise );
% 
% N_data = N_data + length( path_data.mb_0p10_nonoise );
% 
% return

%==================================================================================================
%                                           Everything, No Noise
%==================================================================================================
%
%
%
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

N_data = length( path_data.stoch_bgd ) + N_data;

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

N_data = length( path_data.mb_0p01 ) + N_data;

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

N_data = length( path_data.mb_0p05 ) + N_data;

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

N_data = length( path_data.mb_0p10 ) + N_data;

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

N_data = length( path_data.mb_0p20 ) + N_data;

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

N_data = length( path_data.mb_0p33 ) + N_data;

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

N_data = length( path_data.mb_0p50 ) + N_data;

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

N_data = length( path_data.full_gd ) + N_data;

end

