%

%{


cd /net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi/


restoredefaultpath; 
addpath( genpath( pwd ));


clear; close all; rPIEalpha_sample_vs_probe

%}


%====================================================================================================================================================
% Duchii plots 
%=============

%{

log10_sigma_optim = 3.95;

%========

load /net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/proberPIEalpha_vs_samplerPIEalpha/metrics_mb01_rPIE_alpha_T_vs_alpha_phi.mat
thetitle = 'mb01';

duchii_plots_min_max( metrics.mb01_noise_rmbg, log10_sigma_optim, '-', [ 0.0, 0.0, 0.0 ], thetitle );

%========

% load /net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/proberPIEalpha_vs_samplerPIEalpha/metrics_mb05_rPIE_alpha_T_vs_alpha_phi_NEW.mat
% thetitle = 'mb05';
% 
% duchii_plots_min_max( metrics.mb05_noise_rmbg, log10_sigma_optim, '-', [ 0.0, 0.0, 0.0 ], thetitle );

%========

% load /net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/proberPIEalpha_vs_samplerPIEalpha/metrics_mb10_rPIE_alpha_T_vs_alpha_phi.mat
% thetitle = 'mb10';
% 
% duchii_plots_min_max( metrics.mb10_noise_rmbg, log10_sigma_optim, '-', [ 0.0, 0.0, 0.8 ], thetitle );

%========

% load /net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/proberPIEalpha_vs_samplerPIEalpha/metrics_mb20_rPIE_alpha_T_vs_alpha_phi.mat
% thetitle = 'mb20';
% 
% duchii_plots_min_max( metrics.mb20_noise_rmbg, log10_sigma_optim, '-', [ 0.0, 0.0, 0.8 ], thetitle );

%========

% load /net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/proberPIEalpha_vs_samplerPIEalpha/metrics_mb50_rPIE_alpha_T_vs_alpha_phi.mat
% thetitle = 'mb50';
% 
% duchii_plots_min_max( metrics.mb50_noise_rmbg, log10_sigma_optim, '-', [ 0.8, 0.0, 0.0 ], thetitle );

%========

% load /net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/proberPIEalpha_vs_samplerPIEalpha/metrics_fullGD_rPIE_alpha_T_vs_alpha_phi.mat
% thetitle = 'fullGD';
% 
% duchii_plots_min_max( metrics.fullGD_noise_rmbg, log10_sigma_optim, '-', [ 0.0, 0.8, 0.0 ], thetitle );

%========

% load /net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/proberPIEalpha_vs_samplerPIEalpha/metrics_blockstochGD_rPIE_alpha_T_vs_alpha_phi.mat
% thetitle = 'blockstochGD';
% 
% duchii_plots_min_max( metrics.blockstochGD_noise_rmbg, log10_sigma_optim, '-', [ 0.8, 0.0, 0.8 ], thetitle );





% % duchii_plots_variance( metrics );

return

%}

%====================================================================================================================================================

% rootpath_data = '';

% rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/proberPIEalpha_vs_samplerPIEalpha/fullGD/';
% rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/proberPIEalpha_vs_samplerPIEalpha/blockstochgrad/';
% rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/proberPIEalpha_vs_samplerPIEalpha/mb50/';
% rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/proberPIEalpha_vs_samplerPIEalpha/mb20/';
% rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/proberPIEalpha_vs_samplerPIEalpha/mb10/';
% rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/proberPIEalpha_vs_samplerPIEalpha/mb05/';
% rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/proberPIEalpha_vs_samplerPIEalpha/mb01/';

%========

% rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/pecoCSSIsimulatedmultislicedata/results/';

%========

% rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/poisson_update_bench/31Mar2022/';
rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/poisson_update_bench/02Apr2022/';

%========

% rootpath_data = '';

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
 
    for kk = 1 : length( path_data.( path_data_fields{ jj } ))
        
        %========

%         metrics.( path_data_fields{ jj } ){ kk } = load_and_plot_OLD( path_data.( path_data_fields{ jj } ){ kk }, ...
%                                                                       N_trials.( path_data_fields{ jj } )( kk ),  ...
%                                                                       path_data_fields{ jj },                     ... 
%                                                                       probe_scaling( kk )                           );           
                                                              
        %========
        
        metrics.( path_data_fields{ jj } ){ kk } = load_and_plot( path_data.( path_data_fields{ jj } ){ kk }, ...
                                                                  N_trials.( path_data_fields{ jj } )( kk ),  ...
                                                                  path_data_fields{ jj },                     ... 
                                                                  probe_scaling( kk )                           );           
        
        
    
    end

end

return

%====================================================================================================================================================

function metrics = load_and_plot( path_data, N_trials, path_data_fields, probe_scaling )
    
    metrics = {};
    
    Nspos = 699;

    %========
    
    sim_ptycho2DTPA = cell( N_trials, 1 );
    
    for ii = 1 : N_trials

        sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

    end
    
    Nmetrics_epoch = length( sim_ptycho2DTPA{ 1 }.sol.metrics.meas_gauss_intensity );

    %========

    meas_gauss_intensity = zeros( Nmetrics_epoch, N_trials, 'single' );
    meas_gauss_magnitude = zeros( Nmetrics_epoch, N_trials, 'single' );
    meas_poisson         = zeros( Nmetrics_epoch, N_trials, 'single' );
    
    for ii = 1 : N_trials
        
        meas_gauss_intensity( :, ii ) = sim_ptycho2DTPA{ ii }.sol.metrics.meas_gauss_intensity / probe_scaling;     
        meas_gauss_magnitude( :, ii ) = sim_ptycho2DTPA{ ii }.sol.metrics.meas_gauss_magnitude / probe_scaling;    
        meas_poisson( :, ii )         = sim_ptycho2DTPA{ ii }.sol.metrics.meas_poiss           / probe_scaling;    
   
    end
    
%     std(meas_gauss_intensity, [], 1)
%     std(meas_gauss_magnitude, [], 1)
%     std(meas_poisson, [], 1)
    
    %=============================================================
    % shift/scale all the metrics so that they're between [ 0, 1 ]
    %=============================================================

    I_m = sim_ptycho2DTPA{ 1 }.expt.meas.D .^ 2;
    log_I_m = log( I_m );
    log_I_m( isinf( log_I_m )) = 0;
    
    
    meas_poisson_offset =  ( I_m - I_m .* log_I_m );
    meas_poisson_offset = sum( meas_poisson_offset(:) ) / Nspos;
    
%     figure; 
%     hold on
%     semilogy( 10 * (1 : length( meas_poisson )), meas_poisson, 'Linewidth', 3  ); 
%     semilogy(  10 * (1 : length( meas_poisson )), 0 * meas_poisson + meas_poisson_offset, '--', 'linewidth', 2, 'color', [ 0, 0, 0 ] ); 
%     hold off
%     set( gca, 'fontweight', 'bold' )
%     grid on
%     xlim( [ 0, 5001 ] )


    
    
    
    
%     figure; 
%     
%     subplot(131)
%     semilogy( meas_gauss_intensity, 'Linewidth', 3 ); 
%     set( gca, 'fontweight', 'bold' )
%     grid on
%     
%     subplot(132)
%     semilogy( meas_gauss_magnitude, 'Linewidth', 3 ); 
%     set( gca, 'fontweight', 'bold' )
%     grid on
%     
%     subplot(133)
%     semilogy( meas_poisson, 'Linewidth', 3 ); 
%     set( gca, 'fontweight', 'bold' )
%     grid on
    
    
    
    %%{
    
    meas_poisson = meas_poisson - meas_poisson_offset;
    
% %     meas_gauss_intensity = meas_gauss_intensity - min( meas_gauss_intensity(:) );
% %     meas_gauss_magnitude = meas_gauss_magnitude - min( meas_gauss_magnitude(:) );
% %     meas_poisson         = meas_poisson         - min( meas_poisson(:) );
%     
%     meas_gauss_intensity = meas_gauss_intensity / max( meas_gauss_intensity(:) );
%     meas_gauss_magnitude = meas_gauss_magnitude / max( meas_gauss_magnitude(:) );
%     meas_poisson         = meas_poisson         / max( meas_poisson(:) );

%     figure; 
%     plot( log10( 1 + meas_gauss_intensity )); 
%     ylim([0, 1])
%     
%     figure; 
%     plot( log10( 1 + meas_gauss_magnitude )); 
%     ylim([0, 1])
%     
%     figure; 
%     plot( log10( 1 + meas_poisson )); 
%     ylim([0, 1])
% 



  
%     y_lim = [10^-4, 1];
    
    figure; 
    
    subplot(131)
    semilogy( 10 * (1 : length(meas_gauss_intensity)), meas_gauss_intensity, 'Linewidth', 3 ); 
%     ylim( y_lim )
    set( gca, 'fontweight', 'bold' )
    grid on
    legend( 'Location', 'northeast' ) 
    xlim([ 0, 5000 ])
    title( '$ \frac{1}{N_s} \sum_s \sum_{\mathbf{q}} \vert I^m_{ \mathbf{q}, s } -  I^e_{ \mathbf{q}, s } \vert^2 $', ...
           'FontWeight', 'bold', ...
           'FontSize', 16,       ...
           'Interpreter', 'latex' );
       
      
    subplot(132)
    semilogy( 10 * (1 : length(meas_gauss_intensity)), meas_gauss_magnitude, 'Linewidth', 3 ); 
%     ylim( y_lim )
    set( gca, 'fontweight', 'bold' )
    grid on
    legend( 'Location', 'northeast' ) 
    xlim([ 0, 5000 ])
    title( '$ \frac{1}{N_s} \sum_s \sum_{\mathbf{q}} \big\vert \sqrt{I^m_{ \mathbf{q}, s }} -  \sqrt{I^e_{ \mathbf{q}, s }} \big\vert^2 $', ...
           'FontWeight', 'bold', ...
           'FontSize', 16,       ...
           'Interpreter', 'latex' );
      
    subplot(133)
    semilogy( 10 * (1 : length(meas_gauss_intensity)), meas_poisson, 'Linewidth', 3 ); 
%     ylim( y_lim )
    set( gca, 'fontweight', 'bold' )
    grid on
    legend( 'Location', 'northeast' ) 
    xlim([ 0, 5000 ])
    title( '$ \frac{1}{N_s} \sum_s \sum_{\mathbf{q}} I^e_{ \mathbf{q}, s } - I^m_{ \mathbf{q}, s } log( I^e_{ \mathbf{q}, s } ) + C_0 $', ...
           'FontWeight', 'bold', ...
           'FontSize', 16,       ...
           'Interpreter', 'latex' );
       %}
       
   %================================
   % plot derivatives of the metrics
   %================================
   
   y_lim = [10^-1, 10^11];
   
[~, ii ] = min( meas_gauss_intensity( end, : ), [], 2 );
diff_meas_gauss_intensity = diff( meas_gauss_intensity( :, ii ));

[~, ii ] = min( meas_gauss_magnitude( end, : ), [], 2 );
diff_meas_gauss_magnitude = diff( meas_gauss_magnitude( :, ii ) );

[~, ii ] = min( meas_poisson( end, : ), [], 2 );
diff_meas_poisson         = diff( meas_poisson( :, ii ) );
    


   
% 
% diff_meas_gauss_intensity = diff( meas_gauss_intensity );
% 
% diff_meas_gauss_magnitude = diff( meas_gauss_magnitude );
% 
% diff_meas_poisson         = diff( meas_poisson );
%     



    figure; 
    
    subplot(131)
    semilogy( 10 * (1 : length(diff_meas_gauss_intensity)), abs( diff_meas_gauss_intensity ), 'Linewidth', 2 ); 
    ylim( y_lim )
    set( gca, 'fontweight', 'bold' )
    grid on
    legend( 'Location', 'northeast' ) 
    xlim([ 0, 5000 ])
    xlabel('Epoch')
    title( '$  \frac{d}{d~Epoch}  \Big( \left[ \frac{1}{N_s} \sum_s \sum_{\mathbf{q}} \vert I^m_{ \mathbf{q}, s } -  I^e_{ \mathbf{q}, s } \vert^2 \right] \left( Epoch \right) \Big) $', ...
           'FontWeight', 'bold', ...
           'FontSize', 16,       ...
           'Interpreter', 'latex' );
       
      
    subplot(132)
    semilogy( 10 * (1 : length(diff_meas_gauss_magnitude)),  abs( diff_meas_gauss_magnitude ), 'Linewidth', 2 ); 
    ylim( y_lim )
    set( gca, 'fontweight', 'bold' )
    grid on
    legend( 'Location', 'northeast' ) 
    xlim([ 0, 5000 ])
    xlabel('Epoch')
    title( '$ \frac{d}{d~Epoch}  \Big( \left[ \frac{1}{N_s} \sum_s \sum_{\mathbf{q}} \big\vert \sqrt{I^m_{ \mathbf{q}, s }} -  \sqrt{I^e_{ \mathbf{q}, s }} \big\vert^2 \right] \left( Epoch \right) \Big) $', ...
           'FontWeight', 'bold', ...
           'FontSize', 16,       ...
           'Interpreter', 'latex' );
      
    subplot(133)
    semilogy( 10 * (1 : length(diff_meas_poisson)), abs( diff_meas_poisson ), 'Linewidth', 2 ); 
    ylim( y_lim )
    set( gca, 'fontweight', 'bold' )
    grid on
    legend( 'Location', 'northeast' ) 
    xlim([ 0, 5000 ])
    xlabel('Epoch')
    title( '$  \frac{d}{d~Epoch}  \Big( \left[ \frac{1}{N_s} \sum_s \sum_{\mathbf{q}} I^e_{ \mathbf{q}, s } - I^m_{ \mathbf{q}, s } log( I^e_{ \mathbf{q}, s } ) + C_0 \right] \left( Epoch \right) \Big) $', ...
           'FontWeight', 'bold', ...
           'FontSize', 16,       ...
           'Interpreter', 'latex' );
   
   
   
   
   
   
   
   
   
   
   return
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
       
       
       
       
%     figure; plot( meas_poisson ./ meas_gauss_magnitude )
    
    

    
    
%     plot( ax1, y1, linestyle, 'marker', '.', 'markersize', 12, 'linewidth', 2, 'color', shadingcolor )
% 
% % ccolor = [ 0.0, 0.0, 0.8 ];
% 
% x     = 1 : 40;
% xconf = [ x, x( end : -1 : 1) ]; 
% yconf = [ log10_meas_all_avg_max_early, fliplr( log10_meas_all_avg_min_early ) ];
% % yconf = [ log10_meas_all_avg + 0.5 * log10_meas_all_avg_variance, fliplr( log10_meas_all_avg - 0.5 * log10_meas_all_avg_variance ) ];
% 
% hold on
% 
% p = fill( ax1, xconf, yconf, shadingcolor );
% p.FaceColor = shadingcolor;      
% p.EdgeColor = 'none';           
% p.FaceAlpha = 0.2;
% 
% hold off
% 
% return




% std( meas_gauss_intensity( end, :) )
% std( meas_gauss_magnitude( end, :) )
% std( meas_poisson( end, :) )
%     
    

    
    
    

%     
%     log_mult = 3;
%     
%     meas_gauss_intensity = log10( 1 + 10^log_mult * meas_gauss_intensity );
%     meas_gauss_magnitude = log10( 1 + 10^log_mult * meas_gauss_magnitude );
%     meas_poisson         = log10( 1 + 10^log_mult * meas_poisson );
%     
%     meas_gauss_intensity = meas_gauss_intensity / max( meas_gauss_intensity(:) );
%     meas_gauss_magnitude = meas_gauss_magnitude / max( meas_gauss_magnitude(:) );
%     meas_poisson         = meas_poisson         / max( meas_poisson(:) );
%     
%     
%     figure; 
%     plot( meas_gauss_intensity ); 
% 
%     figure; 
%     plot( meas_gauss_magnitude ); 
% 
%     figure; 
%     plot( meas_poisson ); 
% 
%     return

    %========
% 
%     log10_meas_gauss_intensity_avg = 0;
%     log10_meas_gauss_magnitude_avg = 0;
%     log10_meas_poisson_avg         = 0;
%     
%     for ii = 1 : N_trials
%   
%         log10_meas_gauss_intensity_avg = log10_meas_gauss_intensity_avg + log10( 1 + meas_gauss_intensity( :, ii ));    
%         log10_meas_gauss_magnitude_avg = log10_meas_gauss_magnitude_avg + log10( 1 + meas_gauss_magnitude( :, ii ));    
%         log10_meas_poisson_avg         = log10_meas_poisson_avg         + log10( 1 + meas_poisson( :, ii ));    
%   
%     end
% 
%     log10_meas_gauss_intensity_avg = log10_meas_gauss_intensity_avg / N_trials;
%     log10_meas_gauss_magnitude_avg = log10_meas_gauss_magnitude_avg / N_trials;
%     log10_meas_poisson_avg         = log10_meas_poisson_avg         / N_trials;
% 
%     log10_meas_gauss_intensity_avg = 10 .^ log10_meas_gauss_intensity_avg;
%     log10_meas_gauss_magnitude_avg = 10 .^ log10_meas_gauss_magnitude_avg;
%     log10_meas_poisson_avg         = 10 .^ log10_meas_poisson_avg;
    
    %=================================================
    % create a struct from metrics data for future use
    %=================================================
    
    for ii = 1 : N_trials
        
        metrics.meas_gauss_intensity( :, ii ) = sim_ptycho2DTPA{ ii }.sol.metrics.meas_gauss_intensity;
        metrics.meas_gauss_magnitude( :, ii ) = sim_ptycho2DTPA{ ii }.sol.metrics.meas_gauss_magnitude;
        metrics.meas_poisson( :, ii )         = sim_ptycho2DTPA{ ii }.sol.metrics.meas_poiss;
        
        metrics.timing( ii )          = sim_ptycho2DTPA{ ii }.sol.timings;
        metrics.it( ii )              = sim_ptycho2DTPA{ ii }.sol.it;
        metrics.rPIE_alpha_phi( ii )  = sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_phi;
        metrics.rPIE_alpha_T( ii )    = sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_T;
        metrics.probe_intensity( ii ) = sim_ptycho2DTPA{ ii }.sol.probe.scpm.fro2TOT;
        
        if isfield( sim_ptycho2DTPA{ ii }.sol.spos, 'rand_spos_subset_pct' )
            
            metrics.rand_spos_subset_pct( ii ) = sim_ptycho2DTPA{ ii }.sol.spos.rand_spos_subset_pct;
            
        end
        
    end
    
    %============================================
    % create name strings for the saved the image
    %============================================
    
    if isfield( sim_ptycho2DTPA{ ii }.sol.spos, 'rand_spos_subset_pct' )
        
        name_data = [ path_data_fields, num2str( [ metrics.rPIE_alpha_T( 1 ),           ...
                                                   metrics.rand_spos_subset_pct( 1 )  ], ', rPIE_alpha_T = %0.3e, MBpct = %0.4f') ];
                                               
        save_data = [ path_data_fields, num2str( [ metrics.rPIE_alpha_phi( 1 ),         ...
                                                   metrics.rPIE_alpha_T( 1 ),           ...
                                                   metrics.rand_spos_subset_pct( 1 ) ,  ...
                                                   round( metrics.probe_intensity( 1 )) ], '_rPIEalpha_phi_%0.3e_rPIEalpha_T_%0.3e_MBpct_%0.4f_photons_%0.3e_log10mean'), '.jpg' ];
        
    else
        
        name_data = [ path_data_fields, num2str( [ sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_T ], ', rPIE_alpha_T = %0.3e') ];
                                               
        save_data = [ path_data_fields, num2str( [ sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_T, ...
                                                   round( metrics.probe_intensity( 1 )) ], '_rPIEalpha_phi_%0.3e_rPIEalpha_T_%0.3e_photons_%0.3e_log10mean'), '.jpg' ];
        
    end
  
    %========
    
    metrics.meas_gauss_intensity_max = meas_gauss_intensity_max;
    metrics.meas_gauss_intensity_min = meas_gauss_intensity_min;
    
    metrics.meas_gauss_magnitude_max = meas_gauss_magnitude_max;
    metrics.meas_gauss_magnitude_min = meas_gauss_magnitude_min;
    
    metrics.meas_poisson_max         = meas_poisson_max;
    metrics.meas_poisson_min         = meas_poisson_min;
   
    metrics.meas_gauss_intensity = meas_gauss_intensity;
    metrics.meas_gauss_magnitude = meas_gauss_magnitude;
    metrics.meas_poisson         = meas_poisson;
    
    metrics.meas_gauss_intensity_avg = log10_meas_gauss_intensity_avg;
    metrics.meas_gauss_magnitude_avg = log10_meas_gauss_magnitude_avg;
    metrics.meas_poisson_avg         = log10_meas_poisson_avg;

    metrics.name_data           = name_data;
    metrics.N_trials            = N_trials;
    metrics.path_data           = path_data;
    
    %============================================================
    % create figures of the specified metrics over the trials run
    %============================================================
    
%     skip = 1;
%     
%     y_lim = [ 3, 6 ];
%     
%     h1 = figure();  
%     set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
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
%     plot( metrics.it( 1 ).mtot( 1 : skip : end ), log10_meas_gauss_intensity_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
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
    % ?????
    %==========================================================================

    normalize_metric_start = logical( 0 ); %#ok<LOGL>

    if normalize_metric_start == true
        
%         tmp0 = log10( metrics.meas_all / probe_scaling );
%         tmp1 = log10_meas_gauss_intensity_avg / probe_scaling;
%  
%         tmp2 = log10( metrics.meas_all / probe_scaling );
%         tmp3 = log10_meas_gauss_intensity_avg / probe_scaling;
%         
%         tmp4 = log10( metrics.meas_all / probe_scaling );
%         tmp5 = log10_meas_gauss_intensity_avg / probe_scaling;
        
    else

%         tmp0 = log10( metrics.meas_gauss_intensity );
%         tmp1 = metrics.log10_meas_gauss_intensity_avg;
  
        tmp0 = metrics.meas_gauss_intensity;
        tmp1 = metrics.meas_gauss_intensity_avg;
        
%         tmp2 = log10( metrics.meas_gauss_magnitude );
%         tmp3 = metrics.log10_meas_gauss_magnitude_avg;
        
        tmp2 = metrics.meas_gauss_magnitude;
        tmp3 = metrics.meas_gauss_magnitude_avg;
        
        tmp4 = metrics.meas_poisson;
        tmp5 = metrics.meas_poisson_avg;
        
    end

    tmp6 = metrics.it( 1 ).mtot;
    
    %=========================
    % create the plot, save it
    %=========================
    
    skip = 1;
    
%     y_lim = [ 5.0, 6 ];
%     y_lim = 10.^[ 6.5, 11.0 ];
    y_lim = 10.^[ -7, 0 ];
    
    h1 = figure();  
%     set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
    set( h1, 'Visible', 'on', 'Position',[ 10, 1, 700, 1080 ] )

    
    semilogy( tmp6( 1 : skip : end - 0 ), tmp0( 1 : skip : end - 0, : ), '-', 'linewidth', 2 )
    
    hold on
    semilogy( tmp6( 1 : skip : end - 0 ), tmp1( 1 : skip : end - 0 ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
    hold off
    
    set( gca, 'fontweight', 'bold' )
    xlabel('Epoch')
    ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
    
    title('$ \frac{1}{N_s} \sum_s \sum_{\mathbf{r}}  \big\vert I^m_{ \mathbf{r}, s } -  I^e_{ \mathbf{r}, s }  \big\vert^2 $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
    grid on
    ylim( y_lim )
    legend( 'Location', 'northeast' ) 
    
    export_fig( [ 'intensity_metric_', save_data ], '-r120.0' )
    close all;
    
    %========
    
    skip = 1;
    
%     y_lim = [ 5.0, 6 ];
%     y_lim = 10.^[ 3.5, 6.5 ];

    h1 = figure();  
%     set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
    set( h1, 'Visible', 'on', 'Position',[ 10, 1, 700, 1080 ] )
    
    semilogy( tmp6( 1 : skip : end - 0 ), tmp2( 1 : skip : end - 0, : ), '-', 'linewidth', 2 )
    
    hold on
    semilogy( tmp6( 1 : skip : end - 0 ), tmp3( 1 : skip : end - 0 ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
    hold off
    
    
    set( gca, 'fontweight', 'bold' )
    xlabel('Epoch')
    ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )

    title('$ \frac{1}{N_s} \sum_s \sum_{\mathbf{r}} \big\vert \sqrt{I^m_{ \mathbf{r}, s }} -  \sqrt{I^e_{ \mathbf{r}, s }} \big\vert^2 $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
    grid on
    ylim( y_lim )
    legend( 'Location', 'northeast' ) 
    
    hold off

    export_fig( [ 'magnitude_metric_' , save_data ], '-r120.0' )
    close all;
        
    %========
    
    skip = 1;
    
%     y_lim = [ 5.0, 6 ];
%     y_lim = 10.^[ 0, 7 ];

    h1 = figure();  
%     set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
    set( h1, 'Visible', 'on', 'Position',[ 10, 1, 700, 1080 ] )
    
    
    semilogy( tmp6( 1 : skip : end - 0 ), tmp4( 1 : skip : end - 0, : ), '-', 'linewidth', 2 )
    
    hold on
    semilogy( tmp6( 1 : skip : end - 0 ), tmp5( 1 : skip : end - 0 ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
    hold off
    
    
    
    set( gca, 'fontweight', 'bold' )
    xlabel('Epoch')
    ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )

    title('$ \frac{1}{N_s} \sum_s \sum_{\mathbf{r}} I^e_{ \mathbf{r}, s } - I^m_{ \mathbf{r}, s } log\big( I^e_{ \mathbf{r}, s }  \big) $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
    grid on
    ylim( y_lim )
    legend( 'Location', 'northeast' ) 
    

    fprintf( [ 'saving dataset = ', save_data, '\n' ] );
    
    export_fig( [ 'poisson_metric_', save_data ], '-r120.0' )
    close all;
    
    
    %========================
    % plot variance of trials
    %========================
    
%     stdev = log10( std( meas_all, [], 2 ));  
% %     stdev = std( log10_meas_all, [], 2 );  
%     h1 = figure();  
%     set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
%     
%     plot( tmp2, stdev, '-', 'linewidth', 2 )
%     set( gca, 'fontweight', 'bold' )
%     xlabel('Epoch')
%     ylabel( { [ name_data, ',' ], 'Variance of log10( Cost Function Value )' }, 'Interpreter', 'none' )
%     title('$Var\left( log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]\right)$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
%     grid on
%     
%     
%     save_data = [ path_data_fields, num2str( [ sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_T, ...
%                                                sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_phi, ...
%                                                round( metrics.probe_intensity( 1 )) ], '_rPIE_alpha_phi_%0.3e_rPIE_alpha_T_%0.3e_photons_%0.3e_log10variance'), '.jpg' ];
%     
%     export_fig( save_data, '-r120.0' )
%     close all;
    
    %============================================================================
    % get trial with lowest final cost, and plot corresponding sample/probe modes
    %============================================================================
    
%     [~, ii ] = min( tmp0( end, : ), [], 2 );
%     
%    
%     [ scpm ] = compute_scpm_photonocc( sim_ptycho2DTPA{ ii }.sol.probe.phi );
% 
%     h1 = figure();        
%     set( h1, 'Visible', 'off', 'Position', [ 1, 1, 1920, 1080 ] )
% 
%     for pp = 1 : sim_ptycho2DTPA{ ii }.sol.probe.scpm.N
% 
%         subplot( 1, double( sim_ptycho2DTPA{ ii }.sol.probe.scpm.N ), double( pp ) )   
%         imagescHSV( sim_ptycho2DTPA{ ii }.sol.probe.phi( :, :, pp ) ); 
% 
%         daspect( [ 1, 1, 1 ]); 
% 
%         title( { num2str(  scpm.occ( pp ), 'occupancy = %.4f' ), ...
%                   num2str(  scpm.fro2TOT, 'fro2TOT = %.4e' ) })
%         grid on;
%         set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
% 
%     end
%         
%     save_data = [ path_data_fields,                                                                                                                                   ...
%                   num2str( [ sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_phi, sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_T, round( metrics.probe_intensity( ii )) ], ...
%                   '_rPIEalpha_phi_%0.3e_rPIEalpha_T_%0.3e_photons_%0.3e'),  '_probes.jpg' ];
%     
%     export_fig( save_data, '-r120.0' )
%     close all;
%     
%     %========
% 
%     h1 = figure();        
%     set( h1, 'Visible', 'off', 'Position', [ 1, 1, 1920, 1080 ] )
%     
%     imagescHSV( sim_ptycho2DTPA{ ii }.sol.sample.T )
%     daspect( [ 1, 1, 1 ]); 
%     
%     save_data = [ path_data_fields,                                                                                                                                     ...
%                   num2str( [ sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_phi, sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_T, round( metrics.probe_intensity( ii )) ], ...
%                   '_rPIEalpha_phi_%0.3e_rPIEalpha_T_%0.3e_photons_%0.3e'),    ...
%                   '_sample.jpg' ];
%     
%     export_fig( save_data, '-r120.0' )
%     close all;


end

%====================================================================================================================================================

function metrics = load_and_plot_OLD( path_data, N_trials, path_data_fields, probe_scaling )
    
    metrics = {};

    %========
    
    sim_ptycho2DTPA = cell( N_trials, 1 );
    
    for ii = 1 : N_trials

        sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PetraIII_Eiffel_1p1degree_512x512.mat') ] );
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PetraIII_Eiffel_1p2degree_512x512.mat') ] );
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PetraIII_Eiffel_1p3degree_512x512.mat') ] );
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PetraIII_Eiffel_1p4degree_512x512.mat') ] );

%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PetraIII_Eiffel_1p9degree_512x512.mat') ] );
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PetraIII_Eiffel_1p8degree_512x512.mat') ] );
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PetraIII_Eiffel_1p7degree_512x512.mat') ] );
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PetraIII_Eiffel_1p6degree_512x512.mat') ] );
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PetraIII_Eiffel_1p5degree_512x512.mat') ] );
        
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p1degree_512x512.mat') ] );
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p11degree_512x512.mat') ] );
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p12degree_512x512.mat') ] );  
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p13degree_512x512.mat') ] );
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p14degree_512x512.mat') ] );  
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p15degree_512x512.mat') ] );
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p16degree_512x512.mat') ] );  
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p17degree_512x512.mat') ] );
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p18degree_512x512.mat') ] );  
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p19degree_512x512.mat') ] );  
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p2degree_512x512.mat') ] );
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p3degree_512x512.mat') ] );   
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p4degree_512x512.mat') ] );   
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p5degree_512x512.mat') ] );   
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p6degree_512x512.mat') ] );   
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p7degree_512x512.mat') ] );  
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p8degree_512x512.mat') ] );  
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/PETRAIII_Eiffel_Small_noTi_74ol_1p9degree_512x512.mat') ] );  
        
        
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
    
    %=================================================
    % create a struct from metrics data for future use
    %=================================================
    
    for ii = 1 : N_trials
        
        metrics.meas_all( :, ii )     = sim_ptycho2DTPA{ ii }.sol.metrics.meas_all;
        metrics.timing( ii )          = sim_ptycho2DTPA{ ii }.sol.timings;
        metrics.it( ii )              = sim_ptycho2DTPA{ ii }.sol.it;
        metrics.rPIE_alpha_phi( ii )  = sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_phi;
        metrics.rPIE_alpha_T( ii )    = sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_T;
        metrics.probe_intensity( ii ) = sim_ptycho2DTPA{ ii }.sol.probe.scpm.fro2TOT;
        
        if isfield( sim_ptycho2DTPA{ ii }.sol.spos, 'rand_spos_subset_pct' )
            
            metrics.rand_spos_subset_pct( ii ) = sim_ptycho2DTPA{ ii }.sol.spos.rand_spos_subset_pct;
            
        end
        
    end
    
    %========
    
    
    if isfield( sim_ptycho2DTPA{ ii }.sol.spos, 'rand_spos_subset_pct' )
        
        name_data = [ path_data_fields, num2str( [ metrics.rPIE_alpha_phi( 1 ),         ...
                                                   metrics.rPIE_alpha_T( 1 ),           ...
                                                   metrics.rand_spos_subset_pct( 1 )  ], ', rPIE_alpha_phi = %0.3e, rPIE_alpha_T = %0.3e, MBpct = %0.4f') ];
                                               
        save_data = [ path_data_fields, num2str( [ metrics.rPIE_alpha_phi( 1 ),         ...
                                                   metrics.rPIE_alpha_T( 1 ),           ...
                                                   metrics.rand_spos_subset_pct( 1 ) ,  ...
                                                   round( metrics.probe_intensity( 1 )) ], '_rPIEalpha_phi_%0.3e_rPIEalpha_T_%0.3e_MBpct_%0.4f_photons_%0.3e_log10mean'), '.jpg' ];
        
    else
        
        name_data = [ path_data_fields, num2str( [ sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_phi, ...
                                                   sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_T ], ', rPIE_alpha_phi = %0.3e, rPIE_alpha_T = %0.3e') ];
                                               
        save_data = [ path_data_fields, num2str( [ sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_phi, ...
                                                   sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_T, ...
                                                   round( metrics.probe_intensity( 1 )) ], '_rPIEalpha_phi_%0.3e_rPIEalpha_T_%0.3e_photons_%0.3e_log10mean'), '.jpg' ];
        
    end
    
%     if isfield( sim_ptycho2DTPA{ ii }.sol.spos, 'rand_spos_subset_pct' )
%         
%         name_data = [ path_data_fields, num2str( [ metrics.rPIE_alpha_phi( 1 ),         ...
%                                                    metrics.rPIE_alpha_T( 1 ),           ...
%                                                    metrics.rand_spos_subset_pct( 1 ) ,  ...
%                                                    round( metrics.probe_intensity( 1 )) ], ', rPIE_alpha_phi = %0.3e, rPIE_alpha_T = %0.3e, MBpct = %0.4f, probe photons = %0.3e') ];
%                                                
%         save_data = [ path_data_fields, num2str( [ metrics.rPIE_alpha_phi( 1 ),         ...
%                                                    metrics.rPIE_alpha_T( 1 ),           ...
%                                                    metrics.rand_spos_subset_pct( 1 ) ,  ...
%                                                    round( metrics.probe_intensity( 1 )) ], '_rPIEalpha_phi_%0.3e_rPIEalpha_T_%0.3e_MBpct_%0.4f_photons_%0.3e_log10mean'), '.jpg' ];
%         
%     else
%         
%         name_data = [ path_data_fields, num2str( [ sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_phi, ...
%                                                    sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_T, ...
%                                                    round( metrics.probe_intensity( 1 )) ], ', rPIE_alpha_phi = %0.3e, rPIE_alpha_T = %0.3e, probe photons = %0.3e') ];
%                                                
%         save_data = [ path_data_fields, num2str( [ sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_phi, ...
%                                                    sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_T, ...
%                                                    round( metrics.probe_intensity( 1 )) ], '_rPIEalpha_phi_%0.3e_rPIEalpha_T_%0.3e_photons_%0.3e_log10mean'), '.jpg' ];
%         
%     end
    
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
%     y_lim = [ 3, 6 ];
%     
%     h1 = figure();  
%     set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
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
    
%     y_lim = [ 5.0, 6 ];
    y_lim = [ 3.5, 6.5 ];

    h1 = figure();  
%     set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
    set( h1, 'Visible', 'off', 'Position',[ 10, 1, 700, 1080 ] )
    
    hold on

    for ii = 1 : N_trials

        plot( tmp2( 1 : skip : end ), tmp0( 1 : skip : end, ii ), '-', 'linewidth', 2 )

    end

    plot( tmp2( 1 : skip : end ), tmp1( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
    
%     plot( tmp2( 1 : skip : end ), 4.2 + 0 * tmp1( 1 : skip : end ), '-', 'linewidth', 4, 'color', [ 0.5, 0.5, 0.5 ] )
    
    set( gca, 'fontweight', 'bold' )
    xlabel('Epoch')
    ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
    hold off
    title('$log_{10}\bigg[ \left\langle \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} - \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \right\rangle \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
    grid on
    ylim( y_lim )
    legend( 'Location', 'northeast' ) 
    
    hold off

    fprintf( [ 'saving dataset = ', save_data, '\n' ] );
    
    export_fig( save_data, '-r120.0' );
    close all;
    
    %========================
    % plot variance of trials
    %========================
    
%     stdev = log10( std( meas_all, [], 2 ));  
% %     stdev = std( log10_meas_all, [], 2 );  
%     h1 = figure();  
%     set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
%     
%     plot( tmp2, stdev, '-', 'linewidth', 2 )
%     set( gca, 'fontweight', 'bold' )
%     xlabel('Epoch')
%     ylabel( { [ name_data, ',' ], 'Variance of log10( Cost Function Value )' }, 'Interpreter', 'none' )
%     title('$Var\left( log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]\right)$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
%     grid on
%     
%     
%     save_data = [ path_data_fields, num2str( [ sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_T, ...
%                                                sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_phi, ...
%                                                round( metrics.probe_intensity( 1 )) ], '_rPIE_alpha_phi_%0.3e_rPIE_alpha_T_%0.3e_photons_%0.3e_log10variance'), '.jpg' ];
%     
%     export_fig( save_data, '-r120.0' )
%     close all;
    
    %============================================================================
    % get trial with lowest final cost, and plot corresponding sample/probe modes
    %============================================================================
    
%     [~, ii ] = min( tmp0( end, : ), [], 2 );
%     
%    
%     [ scpm ] = compute_scpm_photonocc( sim_ptycho2DTPA{ ii }.sol.probe.phi );
% 
%     h1 = figure();        
%     set( h1, 'Visible', 'off', 'Position', [ 1, 1, 1920, 1080 ] )
% 
%     for pp = 1 : sim_ptycho2DTPA{ ii }.sol.probe.scpm.N
% 
%         subplot( 1, double( sim_ptycho2DTPA{ ii }.sol.probe.scpm.N ), double( pp ) )   
%         imagescHSV( sim_ptycho2DTPA{ ii }.sol.probe.phi( :, :, pp ) ); 
% 
%         daspect( [ 1, 1, 1 ]); 
% 
%         title( { num2str(  scpm.occ( pp ), 'occupancy = %.4f' ), ...
%                   num2str(  scpm.fro2TOT, 'fro2TOT = %.4e' ) })
%         grid on;
%         set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
% 
%     end
%         
%     save_data = [ path_data_fields,                                                                                                                                   ...
%                   num2str( [ sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_phi, sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_T, round( metrics.probe_intensity( ii )) ], ...
%                   '_rPIEalpha_phi_%0.3e_rPIEalpha_T_%0.3e_photons_%0.3e'),  '_probes.jpg' ];
%     
%     export_fig( save_data, '-r120.0' )
%     close all;
%     
%     %========
% 
%     h1 = figure();        
%     set( h1, 'Visible', 'off', 'Position', [ 1, 1, 1920, 1080 ] )
%     
%     imagescHSV( sim_ptycho2DTPA{ ii }.sol.sample.T )
%     daspect( [ 1, 1, 1 ]); 
%     
%     save_data = [ path_data_fields,                                                                                                                                     ...
%                   num2str( [ sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_phi, sim_ptycho2DTPA{ ii }.sol.rPIE_alpha_T, round( metrics.probe_intensity( ii )) ], ...
%                   '_rPIEalpha_phi_%0.3e_rPIEalpha_T_%0.3e_photons_%0.3e'),    ...
%                   '_sample.jpg' ];
%     
%     export_fig( save_data, '-r120.0' )
%     close all;


end

%====================================================================================================================================================

function [ path_data, N_trials, N_data ] = define_paths( rootpath_data )

N_data = 0;

%================================================

%{

path_data.blockstochGD_noise_rmbg = {};
N_trials.blockstochGD_noise_rmbg  = [];

%========

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-0/independenttrials_27Jan2022_t113930/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-1/independenttrials_02Mar2022_t103812/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-2/independenttrials_03Mar2022_t080005/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-3/independenttrials_28Jan2022_t075512/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-4/independenttrials_04Mar2022_t044619/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-5/independenttrials_05Mar2022_t013001/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-6/independenttrials_29Jan2022_t055218/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-7/independenttrials_05Mar2022_t221253/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-8/independenttrials_06Mar2022_t185706/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-9/independenttrials_30Jan2022_t040931/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

%========

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-0/independenttrials_31Jan2022_t022137/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-1/independenttrials_10Mar2022_t111538/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-2/independenttrials_11Mar2022_t092017/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-3/independenttrials_31Jan2022_t220413/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-4/independenttrials_12Mar2022_t070021/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-5/independenttrials_13Mar2022_t054238/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-6/independenttrials_01Feb2022_t145859/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-7/independenttrials_14Mar2022_t022346/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-8/independenttrials_14Mar2022_t182342/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-9/independenttrials_02Feb2022_t074908/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

%========

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-0/independenttrials_27Jan2022_t114048/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-1/independenttrials_13Mar2022_t044636/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-2/independenttrials_14Mar2022_t020805/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-3/independenttrials_28Jan2022_t080037/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-4/independenttrials_14Mar2022_t183020/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-5/independenttrials_15Mar2022_t105516/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-6/independenttrials_29Jan2022_t060907/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-7/independenttrials_16Mar2022_t031704/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-8/independenttrials_16Mar2022_t193732/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-9/independenttrials_30Jan2022_t043658/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

%========

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-0/independenttrials_31Jan2022_t025800/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-1/independenttrials_07Mar2022_t154222/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-2/independenttrials_08Mar2022_t123335/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-3/independenttrials_31Jan2022_t224258/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-4/independenttrials_09Mar2022_t095431/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-5/independenttrials_10Mar2022_t090458/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-6/independenttrials_01Feb2022_t155109/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-7/independenttrials_11Mar2022_t071121/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-8/independenttrials_12Mar2022_t052751/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

path_data.blockstochGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_blockstochgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-9/independenttrials_02Feb2022_t085909/' ];
N_trials.blockstochGD_noise_rmbg( end + 1 )  = 10;

%========

path_data.blockstochGD_noise_rmbg = transpose( path_data.blockstochGD_noise_rmbg );
N_trials.blockstochGD_noise_rmbg  = transpose( N_trials.blockstochGD_noise_rmbg );

N_data = N_data + length( path_data.blockstochGD_noise_rmbg );

%}

%================================================

%{

path_data.fullGD_noise_rmbg = {};
N_trials.fullGD_noise_rmbg  = [];

%========

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-0/independenttrials_28Jan2022_t091400/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-1/independenttrials_02Mar2022_t115226/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-2/independenttrials_03Mar2022_t032018/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-3/independenttrials_29Jan2022_t005010/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-4/independenttrials_03Mar2022_t184721/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-5/independenttrials_04Mar2022_t101224/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-6/independenttrials_29Jan2022_t162953/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-7/independenttrials_05Mar2022_t013638/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-8/independenttrials_05Mar2022_t170202/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-9/independenttrials_30Jan2022_t081244/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

%========

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-0/independenttrials_04Feb2022_t182612/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-1/independenttrials_09Mar2022_t153018/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-2/independenttrials_10Mar2022_t103605/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-3/independenttrials_06Feb2022_t085140/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-4/independenttrials_11Mar2022_t032523/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-5/independenttrials_11Mar2022_t200936/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-6/independenttrials_07Feb2022_t002631/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-7/independenttrials_12Mar2022_t125927/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-8/independenttrials_13Mar2022_t064241/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-9/independenttrials_07Feb2022_t155252/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

%========


path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-0/independenttrials_08Feb2022_t073002/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-1/independenttrials_09Mar2022_t153510/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-2/independenttrials_10Mar2022_t103939/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-3/independenttrials_08Feb2022_t230718/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-4/independenttrials_11Mar2022_t032758/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-5/independenttrials_11Mar2022_t201122/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-6/independenttrials_09Feb2022_t144354/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-7/independenttrials_12Mar2022_t130054/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-8/independenttrials_13Mar2022_t064316/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-9/independenttrials_10Feb2022_t062020/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

%========

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-0/independenttrials_10Feb2022_t225640/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-1/independenttrials_06Mar2022_t082802/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-2/independenttrials_06Mar2022_t235057/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-3/independenttrials_11Feb2022_t161759/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-4/independenttrials_07Mar2022_t151600/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-5/independenttrials_08Mar2022_t064346/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-6/independenttrials_12Feb2022_t094441/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-7/independenttrials_08Mar2022_t220725/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-8/independenttrials_09Mar2022_t133305/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

path_data.fullGD_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_fullgrad_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-9/independenttrials_13Feb2022_t030938/' ];
N_trials.fullGD_noise_rmbg( end + 1 )  = 10;

%========

path_data.fullGD_noise_rmbg = transpose( path_data.fullGD_noise_rmbg );
N_trials.fullGD_noise_rmbg  = transpose( N_trials.fullGD_noise_rmbg );

N_data = N_data + length( path_data.fullGD_noise_rmbg );

%}

%================================================

%{

% path_data.mb50_noise_rmbg = {};
% N_trials.mb50_noise_rmbg  = [];

%========

% path_data.mb50_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb50_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e0/independenttrials_08Feb2022_t080654/' ];
% N_trials.mb50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb50_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb50_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-3/independenttrials_09Feb2022_t030847/' ];
% N_trials.mb50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb50_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb50_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-6/independenttrials_09Feb2022_t220633/' ];
% N_trials.mb50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb50_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb50_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-9/independenttrials_10Feb2022_t203941/' ];
% N_trials.mb50_noise_rmbg( end + 1 )  = 10;

%========

% path_data.mb50_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb50_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e0/independenttrials_12Feb2022_t082035/' ];
% N_trials.mb50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb50_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb50_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-3/independenttrials_13Feb2022_t200924/' ];
% N_trials.mb50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb50_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb50_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-6/independenttrials_14Feb2022_t180533/' ];
% N_trials.mb50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb50_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb50_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-9/independenttrials_15Feb2022_t070007/' ];
% N_trials.mb50_noise_rmbg( end + 1 )  = 10;

%========

% path_data.mb50_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb50_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e0/independenttrials_07Feb2022_t135436/' ];
% N_trials.mb50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb50_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb50_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-3/independenttrials_08Feb2022_t112148/' ];
% N_trials.mb50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb50_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb50_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-6/independenttrials_09Feb2022_t083527/' ];
% N_trials.mb50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb50_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb50_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-9/independenttrials_10Feb2022_t125802/' ];
% N_trials.mb50_noise_rmbg( end + 1 )  = 10;

%========

% path_data.mb50_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb50_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e0/independenttrials_11Feb2022_t165122/' ];
% N_trials.mb50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb50_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb50_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-3/independenttrials_12Feb2022_t205116/' ];
% N_trials.mb50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb50_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb50_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-6/independenttrials_14Feb2022_t174646/' ];
% N_trials.mb50_noise_rmbg( end + 1 )  = 10;
% 
% path_data.mb50_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb50_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-9/independenttrials_15Feb2022_t051819/' ];
% N_trials.mb50_noise_rmbg( end + 1 )  = 10;

%========

% path_data.mb50_noise_rmbg = transpose( path_data.mb50_noise_rmbg );
% N_trials.mb50_noise_rmbg  = transpose( N_trials.mb50_noise_rmbg );
% 
% N_data = N_data + length( path_data.mb50_noise_rmbg );

%}

%================================================

%{

path_data.mb20_noise_rmbg = {};
N_trials.mb20_noise_rmbg  = [];

%========

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e0/independenttrials_22Jan2022_t143721/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-1/independenttrials_24Mar2022_t141922/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-2/independenttrials_24Mar2022_t231446/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-3/independenttrials_22Jan2022_t143747/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-4/independenttrials_25Mar2022_t080731/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-5/independenttrials_25Mar2022_t170246/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-6/independenttrials_22Jan2022_t143827/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-7/independenttrials_26Mar2022_t015836/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-8/independenttrials_26Mar2022_t105417/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-9/independenttrials_22Jan2022_t143856/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

%========

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e0/independenttrials_23Jan2022_t004333/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-1/independenttrials_26Mar2022_t203214/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-2/independenttrials_27Mar2022_t065201/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-3/independenttrials_23Jan2022_t003310/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-4/independenttrials_27Mar2022_t170624/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-5/independenttrials_28Mar2022_t022519/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-6/independenttrials_23Jan2022_t004452/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-7/independenttrials_28Mar2022_t112706/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-8/independenttrials_28Mar2022_t202957/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-9/independenttrials_23Jan2022_t003424/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

%========

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e0/independenttrials_23Jan2022_t111653/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-1/independenttrials_25Mar2022_t153703/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-2/independenttrials_26Mar2022_t002740/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-3/independenttrials_23Jan2022_t103318/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-4/independenttrials_26Mar2022_t091847/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-5/independenttrials_26Mar2022_t184947/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-6/independenttrials_23Jan2022_t110707/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-7/independenttrials_27Mar2022_t051019/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-8/independenttrials_27Mar2022_t161446/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-9/independenttrials_23Jan2022_t103248/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

%========

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e0/independenttrials_23Jan2022_t214839/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-1/independenttrials_28Mar2022_t020017/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-2/independenttrials_28Mar2022_t112734/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-3/independenttrials_23Jan2022_t203610/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-4/independenttrials_28Mar2022_t205539/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-5/independenttrials_29Mar2022_t061717/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-6/independenttrials_23Jan2022_t213219/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-7/independenttrials_29Mar2022_t150441/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-8/independenttrials_29Mar2022_t235142/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

path_data.mb20_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb20_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-9/independenttrials_23Jan2022_t203555/' ];
N_trials.mb20_noise_rmbg( end + 1 )  = 10;

%========

path_data.mb20_noise_rmbg = transpose( path_data.mb20_noise_rmbg );
N_trials.mb20_noise_rmbg  = transpose( N_trials.mb20_noise_rmbg );

N_data = N_data + length( path_data.mb20_noise_rmbg );

%}

%================================================

%{

path_data.mb10_noise_rmbg = {};
N_trials.mb10_noise_rmbg  = [];

%========

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-0/independenttrials_23Feb2022_t161019/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-1/independenttrials_22Mar2022_t051042/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-2/independenttrials_22Mar2022_t125010/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-3/independenttrials_24Feb2022_t010342/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-4/independenttrials_22Mar2022_t201440/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-5/independenttrials_20Mar2022_t150401/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-6/independenttrials_24Feb2022_t095810/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-7/independenttrials_21Mar2022_t012502/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-8/independenttrials_23Mar2022_t152834/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-9/independenttrials_24Feb2022_t190059/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

%========

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-0/independenttrials_25Feb2022_t040628/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-1/independenttrials_23Mar2022_t230503/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-2/independenttrials_24Mar2022_t063859/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-3/independenttrials_25Feb2022_t131632/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-4/independenttrials_23Mar2022_t033341/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-5/independenttrials_23Mar2022_t105359/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-6/independenttrials_25Feb2022_t222911/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-7/independenttrials_23Mar2022_t181936/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-8/independenttrials_24Mar2022_t014510/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-9/independenttrials_26Feb2022_t073813/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

%========

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-0/independenttrials_26Feb2022_t164036/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-1/independenttrials_20Mar2022_t151500/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-2/independenttrials_21Mar2022_t014221/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-3/independenttrials_27Feb2022_t014632/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-4/independenttrials_23Mar2022_t114019/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-5/independenttrials_23Mar2022_t190437/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-6/independenttrials_27Feb2022_t094638/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-7/independenttrials_24Mar2022_t023431/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-8/independenttrials_24Mar2022_t100602/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-9/independenttrials_27Feb2022_t171338/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

%========

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-0/independenttrials_28Feb2022_t003845/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-1/independenttrials_20Mar2022_t152300/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-2/independenttrials_21Mar2022_t013944/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-3/independenttrials_28Feb2022_t080436/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-4/independenttrials_24Mar2022_t091544/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-5/independenttrials_24Mar2022_t165828/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-6/independenttrials_28Feb2022_t153016/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-7/independenttrials_25Mar2022_t004722/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-8/independenttrials_25Mar2022_t083314/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

path_data.mb10_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb10_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-9/independenttrials_28Feb2022_t225947/' ];
N_trials.mb10_noise_rmbg( end + 1 )  = 10;

%========

path_data.mb10_noise_rmbg = transpose( path_data.mb10_noise_rmbg );
N_trials.mb10_noise_rmbg  = transpose( N_trials.mb10_noise_rmbg );

N_data = N_data + length( path_data.mb10_noise_rmbg );

%}

%================================================

%{

path_data.mb05_noise_rmbg = {};
N_trials.mb05_noise_rmbg  = [];

%========

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e0/independenttrials_24Jan2022_t123102/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-1/independenttrials_17Feb2022_t120722/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-2/independenttrials_18Feb2022_t130421/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-3/independenttrials_24Jan2022_t123142/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-4/independenttrials_19Feb2022_t002331/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-5/independenttrials_19Feb2022_t204236/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-6/independenttrials_14Mar2022_t124507/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-7/independenttrials_20Feb2022_t170056/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-8/independenttrials_21Feb2022_t131213/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e0_rPIE_alpha_T_1e-9/independenttrials_15Mar2022_t104713/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

%========

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e0/independenttrials_24Jan2022_t215351/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-1/independenttrials_21Feb2022_t124445/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-2/independenttrials_21Feb2022_t211353/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-3/independenttrials_24Jan2022_t214820/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-4/independenttrials_22Feb2022_t054131/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-5/independenttrials_22Feb2022_t140319/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-6/independenttrials_14Mar2022_t210100/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-7/independenttrials_22Feb2022_t222001/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-8/independenttrials_23Feb2022_t063013/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-9/independenttrials_15Mar2022_t190230/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

%========

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e0/independenttrials_25Jan2022_t071918/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-1/independenttrials_23Feb2022_t160112/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-2/independenttrials_24Feb2022_t020554/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-3/independenttrials_25Jan2022_t070637/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-4/independenttrials_24Feb2022_t121435/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-5/independenttrials_24Feb2022_t224443/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-6/independenttrials_15Mar2022_t051646/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-7/independenttrials_25Feb2022_t091517/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-8/independenttrials_25Feb2022_t201506/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-9/independenttrials_16Mar2022_t031557/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

%========

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e0/independenttrials_25Jan2022_t164617/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-1/independenttrials_17Feb2022_t121830/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-2/independenttrials_18Feb2022_t131018/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-3/independenttrials_25Jan2022_t162627/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-4/independenttrials_18Feb2022_t213003/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-5/independenttrials_19Feb2022_t080026/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-6/independenttrials_15Mar2022_t133110/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-7/independenttrials_19Feb2022_t183042/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-8/independenttrials_20Feb2022_t050234/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-9_OLD/independenttrials_16Mar2022_t112853/' ];
N_trials.mb05_noise_rmbg( end + 1 )  = 10;

% path_data.mb05_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb05_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-9/independenttrials_03Apr2022_t130725/' ];
% N_trials.mb05_noise_rmbg( end + 1 )  = 10;

%========

path_data.mb05_noise_rmbg = transpose( path_data.mb05_noise_rmbg );
N_trials.mb05_noise_rmbg  = transpose( N_trials.mb05_noise_rmbg );

N_data = N_data + length( path_data.mb05_noise_rmbg );

%}

%================================================

%{

path_data.mb01_noise_rmbg = {};
N_trials.mb01_noise_rmbg  = [];

%========

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-0/independenttrials_18Feb2022_t181309/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-1/independenttrials_20Feb2022_t113828/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-2/independenttrials_22Feb2022_t005622/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-3/independenttrials_23Feb2022_t081819/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-4/independenttrials_24Feb2022_t221049/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-5/independenttrials_16Mar2022_t103201/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-6/independenttrials_17Mar2022_t021326/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-7/independenttrials_17Mar2022_t175304/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-8/independenttrials_18Mar2022_t093244/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-0_rPIE_alpha_T_1e-9/independenttrials_19Mar2022_t011819/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

%========

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-0/independenttrials_26Feb2022_t162109/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-1/independenttrials_27Feb2022_t111556/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-2/independenttrials_28Feb2022_t025940/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-3/independenttrials_28Feb2022_t184057/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-4/independenttrials_01Mar2022_t102503/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-5/independenttrials_21Mar2022_t125021/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-6/independenttrials_22Mar2022_t045243/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-7/independenttrials_21Mar2022_t124706/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-8/independenttrials_22Mar2022_t043934/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-3_rPIE_alpha_T_1e-9/independenttrials_22Mar2022_t201644/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

%========

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-0/independenttrials_26Feb2022_t162624/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-1/independenttrials_27Feb2022_t111234/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-2/independenttrials_28Feb2022_t221801/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-3/independenttrials_01Mar2022_t144044/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-4/independenttrials_21Mar2022_t121037/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-5/independenttrials_17Mar2022_t135250/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-6/independenttrials_18Mar2022_t053836/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-7/independenttrials_18Mar2022_t212817/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-8/independenttrials_19Mar2022_t132505/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-6_rPIE_alpha_T_1e-9/independenttrials_20Mar2022_t051828/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

%========

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-0/independenttrials_19Feb2022_t145554/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-1/independenttrials_21Feb2022_t081519/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-2/independenttrials_22Feb2022_t164104/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-3/independenttrials_24Feb2022_t015549/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-4/independenttrials_25Feb2022_t190500/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-5/independenttrials_18Mar2022_t164902/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-6/independenttrials_19Mar2022_t151959/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-7/independenttrials_20Mar2022_t070950/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-8/independenttrials_20Mar2022_t232319/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

path_data.mb01_noise_rmbg{ end + 1 } = [ rootpath_data, 'probe_mb01_rPIE_alpha_phi_1e-9_rPIE_alpha_T_1e-9/independenttrials_21Mar2022_t154901/' ];
N_trials.mb01_noise_rmbg( end + 1 )  = 10;

%========

path_data.mb01_noise_rmbg = transpose( path_data.mb01_noise_rmbg );
N_trials.mb01_noise_rmbg  = transpose( N_trials.mb01_noise_rmbg );

N_data = N_data + length( path_data.mb01_noise_rmbg );

%}

%====================================================================================================================================================

% path_data.poisson_fixedocc_step5{ 1 } = [ rootpath_data, 'cdi_21Mar_Poisson_FixedOccupancies_Poisson5/independenttrials_23Mar2022_t144843/' ];
% N_trials.poisson_fixedocc_step5( 1 )  = 9;
% 
% path_data.poisson_fixedocc_step10{ 1 } = [ rootpath_data, 'cdi_21Mar_Poisson_FixedOccupancies_Poisson10/independenttrials_23Mar2022_t184612/' ];
% N_trials.poisson_fixedocc_step10( 1 )  = 10;
% 
% path_data.poisson_fixedocc_step25{ 1 } = [ rootpath_data, 'cdi_21Mar_Poisson_FixedOccupancies_Poisson25/independenttrials_23Mar2022_t221537/' ];
% N_trials.poisson_fixedocc_step25( 1 )  = 10;
% 
% path_data.poisson_fixedocc_step50{ 1 } = [ rootpath_data, 'cdi_21Mar_Poisson_FixedOccupancies_Poisson50/independenttrials_24Mar2022_t012936/' ];
% N_trials.poisson_fixedocc_step50( 1 )  = 10;
% 
% path_data.poisson_freeocc_step5{ 1 } = [ rootpath_data, 'cdi_21Mar_Poisson_FreeOccupancies_Poisson5/independenttrials_24Mar2022_t043802/' ];
% N_trials.poisson_freeocc_step5( 1 )  = 10;
% 
% path_data.poisson_freeocc_step10{ 1 } = [ rootpath_data, 'cdi_21Mar_Poisson_FreeOccupancies_Poisson10/independenttrials_24Mar2022_t083445/' ];
% N_trials.poisson_freeocc_step10( 1 )  = 10;
% 
% N_data = 6;

%========

% path_data.gaussian_fixedocc{ 1 } = [ rootpath_data, 'ptycho_mb20_rPIEalphaT_1e-1_gaussian_exp0p5fixedSCPMocc/independenttrials_31Mar2022_t093714/' ];
% N_trials.gaussian_fixedocc( 1 )  = 10;
% 
% path_data.poisson_fixedocc_step1{ 1 } = [ rootpath_data, 'ptycho_mb20_rPIEalphaT_1e-1_poisson1_exp0p5fixedSCPMocc/independenttrials_26Mar2022_t190834/' ];
% N_trials.poisson_fixedocc_step1( 1 )  = 10;
% 
% path_data.poisson_fixedocc_step5{ 1 } = [ rootpath_data, 'ptycho_mb20_rPIEalphaT_1e-1_poisson5_exp0p5fixedSCPMocc/independenttrials_28Mar2022_t191936/' ];
% N_trials.poisson_fixedocc_step5( 1 )  = 10;
% 
% path_data.poisson_fixedocc_step10{ 1 } = [ rootpath_data, 'ptycho_mb20_rPIEalphaT_1e-1_poisson10_exp0p5fixedSCPMocc/independenttrials_28Mar2022_t070834/' ];
% N_trials.poisson_fixedocc_step10( 1 )  = 10;
% 
% path_data.poisson_fixedocc_step25{ 1 } = [ rootpath_data, 'ptycho_mb20_rPIEalphaT_1e-1_poisson25_exp0p5fixedSCPMocc/independenttrials_27Mar2022_t193922/' ];
% N_trials.poisson_fixedocc_step25( 1 )  = 10;
% 
% N_data = 4;

%========

% path_data.gaussian_freeocc_rPIEalphaT_1e0{ 1 } = [ rootpath_data, 'ptycho_mb10_rPIEalphaT_1e-0_gaussian/independenttrials_02Apr2022_t135018/' ];
% N_trials.gaussian_freeocc_rPIEalphaT_1e0( 1 )  = 10;

path_data.poisson_freeocc_step1_rPIEalphaT_1e0{ 1 } = [ rootpath_data, 'ptycho_mb10_rPIEalphaT_1e-0_poisson1/independenttrials_02Apr2022_t135805/' ];
N_trials.poisson_freeocc_step1_rPIEalphaT_1e0( 1 )  = 10;

% path_data.poisson_freeocc_step5_rPIEalphaT_1e0{ 1 } = [ rootpath_data, 'ptycho_mb10_rPIEalphaT_1e-0_poisson5/independenttrials_02Apr2022_t140713/' ];
% N_trials.poisson_freeocc_step5_rPIEalphaT_1e0( 1 )  = 10;

% path_data.poisson_freeocc_step10_rPIEalphaT_1e0{ 1 } = [ rootpath_data, 'ptycho_mb10_rPIEalphaT_1e-0_poisson10/independenttrials_02Apr2022_t141038/' ];
% N_trials.poisson_freeocc_step10_rPIEalphaT_1e0( 1 )  = 10;




% path_data.gaussian_freeocc_rPIEalphaT_1em1{ 1 } = [ rootpath_data, 'ptycho_mb10_rPIEalphaT_1e-1_gaussian/independenttrials_02Apr2022_t213443/' ];
% N_trials.gaussian_freeocc_rPIEalphaT_1em1( 1 )  = 10;

% path_data.poisson_freeocc_step1_rPIEalphaT_1em1{ 1 } = [ rootpath_data, 'ptycho_mb10_rPIEalphaT_1e-1_poisson1/independenttrials_03Apr2022_t024918/' ];
% N_trials.poisson_freeocc_step1_rPIEalphaT_1em1( 1 )  = 10;

% path_data.poisson_freeocc_step5_rPIEalphaT_1em1{ 1 } = [ rootpath_data, 'ptycho_mb10_rPIEalphaT_1e-1_poisson5/independenttrials_02Apr2022_t231648/' ];
% N_trials.poisson_freeocc_step5_rPIEalphaT_1em1( 1 )  = 10;

% path_data.poisson_freeocc_step10_rPIEalphaT_1em1{ 1 } = [ rootpath_data, 'ptycho_mb10_rPIEalphaT_1e-1_poisson10/independenttrials_02Apr2022_t225304/' ];
% N_trials.poisson_freeocc_step10_rPIEalphaT_1em1( 1 )  = 10;



% path_data.gaussian_freeocc_rPIEalphaT_1em2{ 1 } = [ rootpath_data, 'ptycho_mb10_rPIEalphaT_1e-2_gaussian/independenttrials_03Apr2022_t052057/' ];
% N_trials.gaussian_freeocc_rPIEalphaT_1em2( 1 )  = 10;

% path_data.poisson_freeocc_step1_rPIEalphaT_1em2{ 1 } = [ rootpath_data, 'ptycho_mb10_rPIEalphaT_1e-2_poisson1/independenttrials_03Apr2022_t154013/' ];
% N_trials.poisson_freeocc_step1_rPIEalphaT_1em2( 1 )  = 10;

% path_data.poisson_freeocc_step5_rPIEalphaT_1em2{ 1 } = [ rootpath_data, 'ptycho_mb10_rPIEalphaT_1e-2_poisson5/independenttrials_03Apr2022_t082645/' ];
% N_trials.poisson_freeocc_step5_rPIEalphaT_1em2( 1 )  = 10;

% path_data.poisson_freeocc_step10_rPIEalphaT_1em2{ 1 } = [ rootpath_data, 'ptycho_mb10_rPIEalphaT_1e-2_poisson10/independenttrials_03Apr2022_t073122/' ];
% N_trials.poisson_freeocc_step10_rPIEalphaT_1em2( 1 )  = 10;


N_data = 1;

%====================================================================================================================================================

% path_data.eiffel_rocking_1p1 = {};
% N_trials.eiffel_rocking_1p1  = [];
% 
% path_data.eiffel_rocking_1p1{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p1degree_512x512/independenttrials_09Mar2022_t104345/' ];
% N_trials.eiffel_rocking_1p1( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p1 = transpose( path_data.eiffel_rocking_1p1 );
% N_trials.eiffel_rocking_1p1  = transpose( N_trials.eiffel_rocking_1p1 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p1 );

%========


% path_data.eiffel_rocking_1p11 = {};
% N_trials.eiffel_rocking_1p11  = [];
% 
% path_data.eiffel_rocking_1p11{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p11degree_512x512/independenttrials_09Mar2022_t110923/' ];
% N_trials.eiffel_rocking_1p11( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p11 = transpose( path_data.eiffel_rocking_1p11 );
% N_trials.eiffel_rocking_1p11  = transpose( N_trials.eiffel_rocking_1p11 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p11 );

%========

% path_data.eiffel_rocking_1p12 = {};
% N_trials.eiffel_rocking_1p12  = [];
% 
% path_data.eiffel_rocking_1p12{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p12degree_512x512/independenttrials_09Mar2022_t113451/' ];
% N_trials.eiffel_rocking_1p12( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p12 = transpose( path_data.eiffel_rocking_1p12 );
% N_trials.eiffel_rocking_1p12  = transpose( N_trials.eiffel_rocking_1p12 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p12 );

%========

% path_data.eiffel_rocking_1p13 = {};
% N_trials.eiffel_rocking_1p13  = [];
% 
% path_data.eiffel_rocking_1p13{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p13degree_512x512/independenttrials_09Mar2022_t120014/' ];
% N_trials.eiffel_rocking_1p13( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p13 = transpose( path_data.eiffel_rocking_1p13 );
% N_trials.eiffel_rocking_1p13  = transpose( N_trials.eiffel_rocking_1p13 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p13 );

%========


% path_data.eiffel_rocking_1p14 = {};
% N_trials.eiffel_rocking_1p14  = [];
% 
% path_data.eiffel_rocking_1p14{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p14degree_512x512/independenttrials_09Mar2022_t122556/' ];
% N_trials.eiffel_rocking_1p14( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p14 = transpose( path_data.eiffel_rocking_1p14 );
% N_trials.eiffel_rocking_1p14  = transpose( N_trials.eiffel_rocking_1p14 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p14 );

%========

% path_data.eiffel_rocking_1p15 = {};
% N_trials.eiffel_rocking_1p15  = [];
% 
% path_data.eiffel_rocking_1p15{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p15degree_512x512/independenttrials_09Mar2022_t125146/' ];
% N_trials.eiffel_rocking_1p15( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p15 = transpose( path_data.eiffel_rocking_1p15 );
% N_trials.eiffel_rocking_1p15  = transpose( N_trials.eiffel_rocking_1p15 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p15 );

%========

% path_data.eiffel_rocking_1p16 = {};
% N_trials.eiffel_rocking_1p16  = [];
% 
% path_data.eiffel_rocking_1p16{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p16degree_512x512/independenttrials_09Mar2022_t131757/' ];
% N_trials.eiffel_rocking_1p16( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p16 = transpose( path_data.eiffel_rocking_1p16 );
% N_trials.eiffel_rocking_1p16  = transpose( N_trials.eiffel_rocking_1p16 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p16 );

%========

% path_data.eiffel_rocking_1p17 = {};
% N_trials.eiffel_rocking_1p17  = [];
% 
% path_data.eiffel_rocking_1p17{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p17degree_512x512/independenttrials_09Mar2022_t134422/' ];
% N_trials.eiffel_rocking_1p17( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p17 = transpose( path_data.eiffel_rocking_1p17 );
% N_trials.eiffel_rocking_1p17  = transpose( N_trials.eiffel_rocking_1p17 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p17 );

%========

% path_data.eiffel_rocking_1p18 = {};
% N_trials.eiffel_rocking_1p18  = [];
% 
% path_data.eiffel_rocking_1p18{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p18degree_512x512/independenttrials_09Mar2022_t141031/' ];
% N_trials.eiffel_rocking_1p18( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p18 = transpose( path_data.eiffel_rocking_1p18 );
% N_trials.eiffel_rocking_1p18  = transpose( N_trials.eiffel_rocking_1p18 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p18 );

%========

% path_data.eiffel_rocking_1p19 = {};
% N_trials.eiffel_rocking_1p19  = [];
% 
% path_data.eiffel_rocking_1p19{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p19degree_512x512/independenttrials_09Mar2022_t143639/' ];
% N_trials.eiffel_rocking_1p19( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p19 = transpose( path_data.eiffel_rocking_1p19 );
% N_trials.eiffel_rocking_1p19  = transpose( N_trials.eiffel_rocking_1p19 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p19 );

%========

% path_data.eiffel_rocking_1p2 = {};
% N_trials.eiffel_rocking_1p2  = [];
% 
% path_data.eiffel_rocking_1p2{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p2degree_512x512/independenttrials_09Mar2022_t111309/' ];
% N_trials.eiffel_rocking_1p2( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p2 = transpose( path_data.eiffel_rocking_1p2 );
% N_trials.eiffel_rocking_1p2  = transpose( N_trials.eiffel_rocking_1p2 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p2 );

%========

% path_data.eiffel_rocking_1p3 = {};
% N_trials.eiffel_rocking_1p3  = [];
% 
% path_data.eiffel_rocking_1p3{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p3degree_512x512/independenttrials_09Mar2022_t113856/' ];
% N_trials.eiffel_rocking_1p3( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p3 = transpose( path_data.eiffel_rocking_1p3 );
% N_trials.eiffel_rocking_1p3  = transpose( N_trials.eiffel_rocking_1p3 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p3 );

%========

% path_data.eiffel_rocking_1p4 = {};
% N_trials.eiffel_rocking_1p4  = [];
% 
% path_data.eiffel_rocking_1p4{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p4degree_512x512/independenttrials_09Mar2022_t120607/' ];
% N_trials.eiffel_rocking_1p4( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p4 = transpose( path_data.eiffel_rocking_1p4 );
% N_trials.eiffel_rocking_1p4  = transpose( N_trials.eiffel_rocking_1p4 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p4 );

%========

% path_data.eiffel_rocking_1p5 = {};
% N_trials.eiffel_rocking_1p5  = [];
% 
% path_data.eiffel_rocking_1p5{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p5degree_512x512/independenttrials_09Mar2022_t123402/' ];
% N_trials.eiffel_rocking_1p5( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p5 = transpose( path_data.eiffel_rocking_1p5 );
% N_trials.eiffel_rocking_1p5  = transpose( N_trials.eiffel_rocking_1p5 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p5 );

%========

% path_data.eiffel_rocking_1p6 = {};
% N_trials.eiffel_rocking_1p6  = [];
% 
% path_data.eiffel_rocking_1p6{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p6degree_512x512/independenttrials_09Mar2022_t130217/' ];
% N_trials.eiffel_rocking_1p6( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p6 = transpose( path_data.eiffel_rocking_1p6 );
% N_trials.eiffel_rocking_1p6  = transpose( N_trials.eiffel_rocking_1p6 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p6 );

%========

% path_data.eiffel_rocking_1p7 = {};
% N_trials.eiffel_rocking_1p7  = [];
% 
% path_data.eiffel_rocking_1p7{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p7degree_512x512/independenttrials_09Mar2022_t133139/' ];
% N_trials.eiffel_rocking_1p7( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p7 = transpose( path_data.eiffel_rocking_1p7 );
% N_trials.eiffel_rocking_1p7  = transpose( N_trials.eiffel_rocking_1p7 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p7 );

%========

% path_data.eiffel_rocking_1p8 = {};
% N_trials.eiffel_rocking_1p8  = [];
% 
% path_data.eiffel_rocking_1p8{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p8degree_512x512/independenttrials_09Mar2022_t140117/' ];
% N_trials.eiffel_rocking_1p8( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p8 = transpose( path_data.eiffel_rocking_1p8 );
% N_trials.eiffel_rocking_1p8  = transpose( N_trials.eiffel_rocking_1p8 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p8 );

%========

% path_data.eiffel_rocking_1p9 = {};
% N_trials.eiffel_rocking_1p9  = [];
% 
% path_data.eiffel_rocking_1p9{ end + 1 } = [ rootpath_data, '/PETRAIII_Eiffel_Small_noTi_74ol_1p9degree_512x512/independenttrials_09Mar2022_t143228/' ];
% N_trials.eiffel_rocking_1p9( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p9 = transpose( path_data.eiffel_rocking_1p9 );
% N_trials.eiffel_rocking_1p9  = transpose( N_trials.eiffel_rocking_1p9 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p9 );












% path_data.eiffel_rocking_1p1 = {};
% N_trials.eiffel_rocking_1p1  = [];
% 
% path_data.eiffel_rocking_1p1{ end + 1 } = [ rootpath_data, '/PetraIII_Eiffel_1p1degree_512x512/independenttrials_06Mar2022_t141250/' ];
% N_trials.eiffel_rocking_1p1( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p1 = transpose( path_data.eiffel_rocking_1p1 );
% N_trials.eiffel_rocking_1p1 = transpose( N_trials.eiffel_rocking_1p1 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p1 );

%========

% path_data.eiffel_rocking_1p9 = {};
% N_trials.eiffel_rocking_1p9  = [];
% 
% path_data.eiffel_rocking_1p9{ end + 1 } = [ rootpath_data, '/PetraIII_Eiffel_1p9degree_512x512/independenttrials_07Mar2022_t132009/' ];
% N_trials.eiffel_rocking_1p9( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p9 = transpose( path_data.eiffel_rocking_1p9 );
% N_trials.eiffel_rocking_1p9 = transpose( N_trials.eiffel_rocking_1p9 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p9 );

%========

% path_data.eiffel_rocking_1p8 = {};
% N_trials.eiffel_rocking_1p8  = [];
% 
% path_data.eiffel_rocking_1p8{ end + 1 } = [ rootpath_data, '/PetraIII_Eiffel_1p8degree_512x512/independenttrials_07Mar2022_t193610/' ];
% N_trials.eiffel_rocking_1p8( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p8 = transpose( path_data.eiffel_rocking_1p8 );
% N_trials.eiffel_rocking_1p8 = transpose( N_trials.eiffel_rocking_1p8 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p8 );

%========

% path_data.eiffel_rocking_1p7 = {};
% N_trials.eiffel_rocking_1p7  = [];
% 
% path_data.eiffel_rocking_1p7{ end + 1 } = [ rootpath_data, '/PetraIII_Eiffel_1p7degree_512x512/independenttrials_08Mar2022_t012917/' ];
% N_trials.eiffel_rocking_1p7( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p7 = transpose( path_data.eiffel_rocking_1p7 );
% N_trials.eiffel_rocking_1p7 = transpose( N_trials.eiffel_rocking_1p7 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p7 );

%========

% path_data.eiffel_rocking_1p6 = {};
% N_trials.eiffel_rocking_1p6  = [];
% 
% path_data.eiffel_rocking_1p6{ end + 1 } = [ rootpath_data, '/PetraIII_Eiffel_1p6degree_512x512/independenttrials_08Mar2022_t070806/' ];
% N_trials.eiffel_rocking_1p6( end + 1 )  = 6;
% 
% path_data.eiffel_rocking_1p6 = transpose( path_data.eiffel_rocking_1p6 );
% N_trials.eiffel_rocking_1p6 = transpose( N_trials.eiffel_rocking_1p6 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p6 );

%========

% path_data.eiffel_rocking_1p5 = {};
% N_trials.eiffel_rocking_1p5  = [];
% 
% path_data.eiffel_rocking_1p5{ end + 1 } = [ rootpath_data, '/PetraIII_Eiffel_1p5degree_512x512/independenttrials_07Mar2022_t075928/' ];
% N_trials.eiffel_rocking_1p5( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p5 = transpose( path_data.eiffel_rocking_1p5 );
% N_trials.eiffel_rocking_1p5 = transpose( N_trials.eiffel_rocking_1p5 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p5 );

%========

% path_data.eiffel_rocking_1p2 = {};
% N_trials.eiffel_rocking_1p2  = [];
% 
% path_data.eiffel_rocking_1p2{ end + 1 } = [ rootpath_data, '/PetraIII_Eiffel_1p2degree_512x512/independenttrials_06Mar2022_t181224/' ];
% N_trials.eiffel_rocking_1p2( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p2 = transpose( path_data.eiffel_rocking_1p2 );
% N_trials.eiffel_rocking_1p2 = transpose( N_trials.eiffel_rocking_1p2 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p2 );

%========

% path_data.eiffel_rocking_1p3 = {};
% N_trials.eiffel_rocking_1p3  = [];
% 
% path_data.eiffel_rocking_1p3{ end + 1 } = [ rootpath_data, '/PetraIII_Eiffel_1p3degree_512x512/independenttrials_06Mar2022_t223325/' ];
% N_trials.eiffel_rocking_1p3( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p3 = transpose( path_data.eiffel_rocking_1p3 );
% N_trials.eiffel_rocking_1p3 = transpose( N_trials.eiffel_rocking_1p3 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p3 );

%========

% path_data.eiffel_rocking_1p4 = {};
% N_trials.eiffel_rocking_1p4  = [];
% 
% path_data.eiffel_rocking_1p4{ end + 1 } = [ rootpath_data, '/PetraIII_Eiffel_1p4degree_512x512/independenttrials_07Mar2022_t030055/' ];
% N_trials.eiffel_rocking_1p4( end + 1 )  = 10;
% 
% path_data.eiffel_rocking_1p4 = transpose( path_data.eiffel_rocking_1p4 );
% N_trials.eiffel_rocking_1p4 = transpose( N_trials.eiffel_rocking_1p4 );
% 
% N_data = N_data + length( path_data.eiffel_rocking_1p4 );


end

%====================================================================================================================================================

function duchii_plots_min_max( metrics_noise_rmbg, log10_sigma_optim, linestyle, shadingcolor, thetitle )

% log10_sigma_optim = 4.1;

% success_fraction = zeros( 40, 1 );

kk = 1;

Ndata = length( metrics_noise_rmbg );

for ii = 1 : Ndata
        
    %========
    
    rPIE_alpha_T( kk )      = metrics_noise_rmbg{ ii }.rPIE_alpha_T( 1 );   
    rPIE_alpha_phi( kk )    = metrics_noise_rmbg{ ii }.rPIE_alpha_phi( 1 );    
    
    log10_meas_all_avg_end( kk )     = metrics_noise_rmbg{ ii }.log10_meas_all_avg( end );
    log10_meas_all_avg_max_end( kk ) = metrics_noise_rmbg{ ii }.log10_meas_all_max( end );
    log10_meas_all_avg_min_end( kk ) = metrics_noise_rmbg{ ii }.log10_meas_all_min( end );
    log10_meas_all_avg_variance_end( kk ) = std( metrics_noise_rmbg{ ii }.log10_meas_all( end, : ) );

    Nearly_epoch = 5;
    log10_meas_all_avg_early( kk )     = metrics_noise_rmbg{ ii }.log10_meas_all_avg( Nearly_epoch );
    log10_meas_all_avg_max_early( kk ) = metrics_noise_rmbg{ ii }.log10_meas_all_max( Nearly_epoch );
    log10_meas_all_avg_min_early( kk ) = metrics_noise_rmbg{ ii }.log10_meas_all_min( Nearly_epoch );
    log10_meas_all_avg_variance_early( kk ) = std( metrics_noise_rmbg{ ii }.log10_meas_all( Nearly_epoch, : ) );
    
    %========
    
%     epoch_indx = find( metrics_noise_rmbg{ ii }.log10_meas_all_avg < log10_sigma_optim, 1 );
%     
%     tmp0 = metrics_noise_rmbg{ ii }.it( 1 ).mtot( epoch_indx );
% 
%     if isempty( tmp0 )
% 
%         epoch_sigopt_mean( kk ) = metrics_noise_rmbg{ ii }.it( 1 ).mtot( end );
%         
%         
%     else
%         
%         epoch_sigopt_mean( kk ) = metrics_noise_rmbg{ ii }.it( 1 ).mtot( epoch_indx ); 
%     
%     end

    %========
    
%     sf = 0;
%     
%     for cc = 1 : 5
% 
%         II = find( metrics_noise_rmbg{ ii }.log10_meas_all( :, cc ) < log10_sigma_optim, 1 );
% 
%         if isempty( II )
%             
% %             success_fraction( kk ) = sf;
% 
%             tmp0( cc ) = metrics_noise_rmbg{ ii }.it( 1 ).mtot( end );
% 
%         else
% 
%             sf = sf + 1;        
%             success_fraction( kk ) = sf;
%             
%             tmp0( cc ) = metrics_noise_rmbg{ ii }.it( 1 ).mtot( II );
% 
%         end
% 
%     end
% 
%     min_epoch_sigopt_10trials( kk ) = min( tmp0 );
%     max_epoch_sigopt_10trials( kk ) = max( tmp0 );
    
    %========
    
%     figure; 
%     plot(metrics_noise_rmbg{ ii }.log10_meas_all); 
%     hold on; 
%     plot(metrics_noise_rmbg{ ii }.log10_meas_all_avg, '--k', 'linewidth', 2 ); 
%     plot( log10_sigma_optim + 0 * metrics_noise_rmbg{ ii }.log10_meas_all, '-k', 'linewidth', 2 ); 
%     hold off;
%     title( num2str( [ rPIE_alpha_T( kk ), rPIE_alpha_phi( kk ) , success_fraction( kk ) ], 'rPIE_alpha_T = %.3e, rPIE_alpha_phi = %.3e, successful = %d' ),'Interpreter','none' )
%     ylim([ 3.5, 6.5 ])
%            
       
    
    %========
    
    kk = kk + 1;

end

%========

for ii = 1 : Ndata
    
    XTickLabel_T{ ii }   = num2str( rPIE_alpha_T( ii ),   '%.1e' );
    XTickLabel_phi{ ii } = num2str( rPIE_alpha_phi( ii ), '%.1e' );
    
end

%========

% x1 = rPIE_alpha_T;
% x2 = rPIE_alpha_phi;
% y1 = log10_meas_all_avg_end;

figure( 1 );
set( gcf, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
% set( gcf, 'Visible', 'on', 'Position',[ 10, 1, 700, 1080 ] )

t = tiledlayout(1,1);
ax1 = axes(t);

%========

slc_alphaphi_1em0 = 1 : 10;
slc_alphaphi_1em3 = 11 : 20;
slc_alphaphi_1em6 = 21 : 30;
slc_alphaphi_1em9 = 31 : 40;

% slc_alphaphi_1em0 = 1 : 10;
% slc_alphaphi_1em3 = 11 : 19;
% slc_alphaphi_1em6 = 20 : 29;
% slc_alphaphi_1em9 = 30 : 39;

% slc_alphaphi_1em0 = 1 : 4;
% slc_alphaphi_1em3 = 5 : 8;
% slc_alphaphi_1em6 = 9 : 12;
% slc_alphaphi_1em9 = 13 : 16;

hold on

plot( ax1, slc_alphaphi_1em0, log10_meas_all_avg_end( slc_alphaphi_1em0 ), linestyle, 'marker', '.', 'markersize', 12, 'linewidth', 2, 'color', [ 0.8, 0.0, 0.0 ] )
plot( ax1, slc_alphaphi_1em3, log10_meas_all_avg_end( slc_alphaphi_1em3 ), linestyle, 'marker', '.', 'markersize', 12, 'linewidth', 2, 'color', [ 0.0, 0.0, 0.8 ] )
plot( ax1, slc_alphaphi_1em6, log10_meas_all_avg_end( slc_alphaphi_1em6 ), linestyle, 'marker', '.', 'markersize', 12, 'linewidth', 2, 'color', [ 0.0, 0.8, 0.0 ] )
plot( ax1, slc_alphaphi_1em9, log10_meas_all_avg_end( slc_alphaphi_1em9 ), linestyle, 'marker', '.', 'markersize', 12, 'linewidth', 2, 'color', [ 0.8, 0.0, 0.8 ] )

plot( ax1, slc_alphaphi_1em0, log10_meas_all_avg_early( slc_alphaphi_1em0 ), '--', 'marker', '.', 'markersize', 12, 'linewidth', 2, 'color', [ 0.8, 0.0, 0.0 ] )
plot( ax1, slc_alphaphi_1em3, log10_meas_all_avg_early( slc_alphaphi_1em3 ), '--', 'marker', '.', 'markersize', 12, 'linewidth', 2, 'color', [ 0.0, 0.0, 0.8 ] )
plot( ax1, slc_alphaphi_1em6, log10_meas_all_avg_early( slc_alphaphi_1em6 ), '--', 'marker', '.', 'markersize', 12, 'linewidth', 2, 'color', [ 0.0, 0.8, 0.0 ] )
plot( ax1, slc_alphaphi_1em9, log10_meas_all_avg_early( slc_alphaphi_1em9 ), '--', 'marker', '.', 'markersize', 12, 'linewidth', 2, 'color', [ 0.8, 0.0, 0.8 ] )

%========

x     = slc_alphaphi_1em0;
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ log10_meas_all_avg_max_end( x ), fliplr( log10_meas_all_avg_min_end( x ) ) ];
p = fill( ax1, xconf, yconf, [ 0.8, 0.0, 0.0 ] );
p.FaceColor = [ 0.8, 0.0, 0.0 ];    
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

x     = slc_alphaphi_1em3;
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ log10_meas_all_avg_max_end( x ), fliplr( log10_meas_all_avg_min_end( x ) ) ];
p = fill( ax1, xconf, yconf, [ 0.0, 0.0, 0.8 ] );
p.FaceColor = [ 0.0, 0.0, 0.8 ];    
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

x     = slc_alphaphi_1em6;
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ log10_meas_all_avg_max_end( x ), fliplr( log10_meas_all_avg_min_end( x ) ) ];
p = fill( ax1, xconf, yconf, [ 0.0, 0.8, 0.0 ] );
p.FaceColor = [ 0.0, 0.8, 0.0 ];    
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

x     = slc_alphaphi_1em9;
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ log10_meas_all_avg_max_end( x ), fliplr( log10_meas_all_avg_min_end( x ) ) ];
p = fill( ax1, xconf, yconf, [ 0.8, 0.0, 0.8 ] );
p.FaceColor = [ 0.8, 0.0, 0.8 ];    
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

%========

x     = slc_alphaphi_1em0;
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ log10_meas_all_avg_max_early( x ), fliplr( log10_meas_all_avg_min_early( x ) ) ];
p = fill( ax1, xconf, yconf, [ 0.4, 0.0, 0.0 ] );
p.FaceColor = [ 0.4, 0.0, 0.0 ];    
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

x     = slc_alphaphi_1em3;
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ log10_meas_all_avg_max_early( x ), fliplr( log10_meas_all_avg_min_early( x ) ) ];
p = fill( ax1, xconf, yconf, [ 0.0, 0.0, 0.4 ] );
p.FaceColor = [ 0.0, 0.0, 0.4 ];    
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

x     = slc_alphaphi_1em6;
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ log10_meas_all_avg_max_early( x ), fliplr( log10_meas_all_avg_min_early( x ) ) ];
p = fill( ax1, xconf, yconf, [ 0.0, 0.4, 0.0 ] );
p.FaceColor = [ 0.0, 0.4, 0.0 ];    
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

x     = slc_alphaphi_1em9;
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ log10_meas_all_avg_max_early( x ), fliplr( log10_meas_all_avg_min_early( x ) ) ];
p = fill( ax1, xconf, yconf, [ 0.4, 0.0, 0.4 ] );
p.FaceColor = [ 0.4, 0.0, 0.4 ];    
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

hold off

%========

ax1.XTick = 1 : Ndata;
% ax1.XLabel.String = 'rPIE\_alpha\_T';
ax1.XLabel.String = '\alpha_T';
ax1.XTickLabel = XTickLabel_T;
ax1.XLim = [ 1, Ndata ];
ax1.YLim = [ 3.5, 6.5 ];
ax1.FontWeight = 'bold';

grid on
title( thetitle, '500th and 5000th Epoch' )

export_fig( [ thetitle, '_500_and_5000.png'], '-r120.0' )
print( [ thetitle, '_500_and_5000.svg' ], '-dsvg' )

close all;

%========

% x1 = rPIE_alpha_T;
% x2 = rPIE_alpha_phi;
% y1 = log10_meas_all_avg_early;
% 
% figure( 1 );
% set( gcf, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
% % set( gcf, 'Visible', 'on', 'Position',[ 10, 1, 700, 1080 ] )
% 
% t = tiledlayout(1,1);
% ax1 = axes(t);
% 
% 
% 
% hold on
% 
% x     = 1 : 10;
% xconf = [ x, x( end : -1 : 1) ]; 
% yconf = [ log10_meas_all_avg_max_early( x ), fliplr( log10_meas_all_avg_min_early( x ) ) ];
% 
% p = fill( ax1, xconf, yconf, [ 0.8, 0.0, 0.0 ] );
% p.FaceColor = [ 0.8, 0.0, 0.0 ];    
% p.EdgeColor = 'none';           
% p.FaceAlpha = 0.2;
% 
% x     = 11 : 20;
% xconf = [ x, x( end : -1 : 1) ]; 
% yconf = [ log10_meas_all_avg_max_early( x ), fliplr( log10_meas_all_avg_min_early( x ) ) ];
% 
% p = fill( ax1, xconf, yconf, [ 0.0, 0.0, 0.8 ] );
% p.FaceColor = [ 0.0, 0.0, 0.8 ];    
% p.EdgeColor = 'none';           
% p.FaceAlpha = 0.2;
% 
% x     = 21 : 30;
% xconf = [ x, x( end : -1 : 1) ]; 
% yconf = [ log10_meas_all_avg_max_early( x ), fliplr( log10_meas_all_avg_min_early( x ) ) ];
% 
% p = fill( ax1, xconf, yconf, [ 0.0, 0.8, 0.0 ] );
% p.FaceColor = [ 0.0, 0.8, 0.0 ];    
% p.EdgeColor = 'none';           
% p.FaceAlpha = 0.2;
% 
% x     = 31 : 40;
% xconf = [ x, x( end : -1 : 1) ]; 
% yconf = [ log10_meas_all_avg_max_early( x ), fliplr( log10_meas_all_avg_min_early( x ) ) ];
% 
% p = fill( ax1, xconf, yconf, [ 0.0, 0.0, 0.0 ] );
% p.FaceColor = [ 0.0, 0.0, 0.0 ];    
% p.EdgeColor = 'none';           
% p.FaceAlpha = 0.2;
% 
% hold off
% 
% ax1.XTick = 1 : 40;
% % ax1.XLabel.String = 'rPIE\_alpha\_T';
% ax1.XLabel.String = '\alpha_T';
% ax1.XTickLabel = XTickLabel_T;
% ax1.XLim = [ 1, 40 ];
% ax1.YLim = [ 3.5, 6.5 ];
% ax1.FontWeight = 'bold';
% 
% grid on
% title( thetitle, '500th Epoch' )
% 
% export_fig( [ thetitle, '_500.png'], '-r120.0' )
% print( [ thetitle, '_500.svg' ], '-dsvg' )
% 
% close all;





















return















































% x1 = rPIE_alpha_T;
% x2 = rPIE_alpha_phi;
% y1 = log10_meas_all_avg_end;
% 
% 
% figure( 1 );
% set( gcf, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
% % set( gcf, 'Visible', 'on', 'Position',[ 10, 1, 700, 1080 ] )
% 
% t = tiledlayout(1,1);
% ax1 = axes(t);
% 
% plot( ax1, y1, linestyle, 'marker', '.', 'markersize', 12, 'linewidth', 2, 'color', shadingcolor )
% 
% % ccolor = [ 0.0, 0.0, 0.8 ];
% 
% x     = 1 : 40;
% xconf = [ x, x( end : -1 : 1) ]; 
% yconf = [ log10_meas_all_avg_max_end, fliplr( log10_meas_all_avg_min_end ) ];
% % yconf = [ log10_meas_all_avg + 0.5 * log10_meas_all_avg_variance, fliplr( log10_meas_all_avg - 0.5 * log10_meas_all_avg_variance ) ];
% 
% hold on
% 
% p = fill( ax1, xconf, yconf, shadingcolor );
% p.FaceColor = shadingcolor;      
% p.EdgeColor = 'none';           
% p.FaceAlpha = 0.2;
% 
% hold off
% 
% ax1.XTick = 1 : 40;
% % ax1.XLabel.String = 'rPIE\_alpha\_T';
% ax1.XLabel.String = '\alpha_T';
% ax1.XTickLabel = XTickLabel_T;
% ax1.XLim = [ 1, 40 ];
% ax1.YLim = [ 3.5, 6.5 ];
% ax1.FontWeight = 'bold';
% 
% 
% ax2 = axes(t);
% plot( ax2, y1, linestyle, 'marker', '.', 'markersize', 12, 'linewidth', 2, 'color', shadingcolor )
% 
% hold on
% 
% p = fill( ax2, xconf, yconf, shadingcolor );
% p.FaceColor = shadingcolor;      
% p.EdgeColor = 'none';           
% p.FaceAlpha = 0.2;
% 
% hold off
% 
% ax2.XAxisLocation = 'top';
% ax2.XTick = 1 : 40;
% % ax2.XLabel.String = 'rPIE\_alpha\_phi';
% ax2.XLabel.String = '\alpha_{\phi}';
% ax2.XTickLabel = XTickLabel_phi;
% ax2.XLim = [ 1 , 40 ];
% ax2.YLim = [ 3.5, 6.5 ];
% ax2.FontWeight = 'bold';
% 
% 
% grid on
% title( thetitle, '5000th Epoch' )
% 
% 
% export_fig( [ thetitle, '_5000.png'], '-r120.0' )
% 
% print( [ thetitle, '_5000.svg' ], '-dsvg' )
% 
% close all;

%========


x1 = rPIE_alpha_T;
x2 = rPIE_alpha_phi;
y1 = log10_meas_all_avg_early;


figure( 2 );
set( gcf, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
% set( gcf, 'Visible', 'on', 'Position',[ 10, 1, 700, 1080 ] )

t = tiledlayout(1,1);
ax1 = axes(t);

plot( ax1, y1, linestyle, 'marker', '.', 'markersize', 12, 'linewidth', 2, 'color', shadingcolor )

% ccolor = [ 0.0, 0.0, 0.8 ];

x     = 1 : 40;
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ log10_meas_all_avg_max_early, fliplr( log10_meas_all_avg_min_early ) ];
% yconf = [ log10_meas_all_avg + 0.5 * log10_meas_all_avg_variance, fliplr( log10_meas_all_avg - 0.5 * log10_meas_all_avg_variance ) ];

hold on

p = fill( ax1, xconf, yconf, shadingcolor );
p.FaceColor = shadingcolor;      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

hold off

ax1.XTick = 1 : 40;
% ax1.XLabel.String = 'rPIE\_alpha\_T';
ax1.XLabel.String = '\alpha_T';
ax1.XTickLabel = XTickLabel_T;
ax1.XLim = [ 1 , 40 ];
ax1.YLim = [ 3.5, 6.5 ];
ax1.FontWeight = 'bold';


ax2 = axes(t);
plot( ax2, y1, linestyle, 'marker', '.', 'markersize', 12, 'linewidth', 2, 'color', shadingcolor )

hold on

p = fill( ax2, xconf, yconf, shadingcolor );
p.FaceColor = shadingcolor;      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

hold off

ax2.XAxisLocation = 'top';
ax2.XTick = 1 : 40;
% ax2.XLabel.String = 'rPIE\_alpha\_phi';
ax2.XLabel.String = '\alpha_{\phi}';
ax2.XTickLabel = XTickLabel_phi;
ax2.XLim = [ 1, 40 ];
ax2.YLim = [ 3.5, 6.5 ];
ax2.FontWeight = 'bold';

grid on
title( thetitle, '500th Epoch' )

export_fig( [ thetitle, '_500.png'], '-r120.0' )

print( [ thetitle, '_500.svg' ], '-dsvg' )

close all;


%========










return



% 
% 
% ccolor = [ 0.0, 0.0, 0.8 ];
% 
% x     = 1 : 16;
% xconf = [ x, x( end : -1 : 1) ]; 
% yconf = [ log10_meas_all_avg_max, fliplr( log10_meas_all_avg_min ) ];
% 
% 
% ax3 = axes(t);
% 
% hold on
% 
% p = fill( ax3, xconf, yconf, ccolor );
% p.FaceColor = ccolor;      
% p.EdgeColor = 'none';           
% p.FaceAlpha = 0.2;
% 
% hold off
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% tmp95 = rPIE_alpha_T;
% tmp96 = rPIE_alpha_phi;
% tmp97 = log10_meas_all_avg; 
% 
% tmp95 = transpose( reshape( tmp95, [ 4, 4 ] )); rPIE_alpha_T       = tmp95(:); 
% tmp96 = transpose( reshape( tmp96, [ 4, 4 ] )); rPIE_alpha_phi     = tmp96(:); 
% tmp97 = transpose( reshape( tmp97, [ 4, 4 ] )); log10_meas_all_avg = tmp97(:); 
% 
% 
% figure
% 
% x1 = rPIE_alpha_T;
% x2 = rPIE_alpha_phi;
% y1 = log10_meas_all_avg;
% 
% 
% t = tiledlayout(1,1);
% 
% ax1 = axes(t);
% plot( ax1, y1, '-r.' )
% ax1.XTick = 1 : 16;
% ax1.XLabel.String = 'rPIE\_alpha\_T';
% ax1.XTickLabel = XTickLabel_T;
% 
% 
% ax2 = axes(t);
% plot( ax2, y1, '-r.' )
% ax2.XAxisLocation = 'top';
% ax2.XTick = 1 : 16;
% ax2.XLabel.String = 'rPIE\_alpha\_phi';
% ax2.XTickLabel = XTickLabel_phi;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % ax1.Color = 'none';
% ax1.XColor = [ 0.8, 0.0, 0.0 ];
% ax1.YColor = [ 0.8, 0.0, 0.0 ];
% 
% ax2 = axes(t);
% plot(ax2,x2,y2,'-k')
% 
% ax2.XAxisLocation = 'top';
% ax2.YAxisLocation = 'right';
% ax2.YLabel.String = 'AAA';
% ax2.XLabel.String = 'BBB';
% 
% 
% ax2.Color = 'none';
% % ax2.Color = [ 0.0, 0.7, 0.0 ];
% ax2.XColor = [ 0.0, 0.7, 0.0 ];
% ax2.YColor = [ 0.0, 0.7, 0.0 ];
% 
% ax2.XLabel.Color = [ 0.0, 0.0, 0.5 ];
% 
% ax1.Box = 'off';
% ax2.Box = 'off';
% 










%========

h1 = figure;  
set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )

%========

data_slice = 1 : 4;
ccolor = [ 0.0, 0.0, 0.8 ];

x     = log10( rPIE_alpha_T( data_slice ));
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ max_epoch_sigopt_10trials( data_slice ), fliplr( min_epoch_sigopt_10trials( data_slice ) ) ];

hold on

p = fill( xconf, yconf, ccolor );
p.FaceColor = ccolor;      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

plot( log10( rPIE_alpha_T( data_slice ) ), epoch_sigopt_mean( data_slice ), '-', 'Marker', '.', 'linewidth', 2, 'markersize', 30, 'color', ccolor )

hold off

grid on
set( gca, 'fontweight', 'bold' )
ylabel('Epoch')
xlabel('log10( rPIE_alpha_T )','Interpreter','none')

% title( num2str( [ log10_sigma_optim, rPIE_alpha_phi( data_slice( end )) ], 'log10_sigma_optim  = %.3f, rPIE_alpha_phi = %.1e'),'Interpreter','none')
title( num2str( log10_sigma_optim, 'log10_sigma_optim  = %.3f'),'Interpreter','none')
% ylim([ 1, 5000 ])

%========

data_slice = 5 : 8;
ccolor = [ 0.0, 0.8, 0.0 ];

x     = log10( rPIE_alpha_T( data_slice ));
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ max_epoch_sigopt_10trials( data_slice ), fliplr( min_epoch_sigopt_10trials( data_slice ) ) ];

hold on

p = fill( xconf, yconf, ccolor );
p.FaceColor = ccolor;      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

plot( log10( rPIE_alpha_T( data_slice ) ), epoch_sigopt_mean( data_slice ), '-', 'Marker', '.', 'linewidth', 2, 'markersize', 30, 'color', ccolor )

hold off

grid on
set( gca, 'fontweight', 'bold' )
ylabel('Epoch')
xlabel('log10( rPIE_alpha_T )','Interpreter','none')
title( num2str( log10_sigma_optim, 'log10_sigma_optim  = %.3f'),'Interpreter','none')
% title( num2str( [ log10_sigma_optim, rPIE_alpha_phi( data_slice( end )) ], 'log10_sigma_optim  = %.3f, rPIE_alpha_phi = %.1e'),'Interpreter','none')
% ylim([ 1, 5000 ])

%========

data_slice = 9 : 12;
ccolor = [ 0.8, 0.0, 0.0 ];

x     = log10( rPIE_alpha_T( data_slice ));
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ max_epoch_sigopt_10trials( data_slice ), fliplr( min_epoch_sigopt_10trials( data_slice ) ) ];

hold on

p = fill( xconf, yconf, ccolor );
p.FaceColor = ccolor;      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

plot( log10( rPIE_alpha_T( data_slice ) ), epoch_sigopt_mean( data_slice ), '-', 'Marker', '.', 'linewidth', 2, 'markersize', 30, 'color', ccolor )

hold off

grid on
set( gca, 'fontweight', 'bold' )
ylabel('Epoch')
xlabel('log10( rPIE_alpha_T )','Interpreter','none')
title( num2str( log10_sigma_optim, 'log10_sigma_optim  = %.3f'),'Interpreter','none')
% title( num2str( [ log10_sigma_optim, rPIE_alpha_phi( data_slice( end )) ], 'log10_sigma_optim  = %.3f, rPIE_alpha_phi = %.1e'),'Interpreter','none')
% ylim([ 1, 6000 ])

%========

data_slice = 13 : 16;
ccolor = [ 0.8, 0.4, 0.8 ];

x     = log10( rPIE_alpha_T( data_slice ));
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ max_epoch_sigopt_10trials( data_slice ), fliplr( min_epoch_sigopt_10trials( data_slice ) ) ];

hold on

p = fill( xconf, yconf, ccolor );
p.FaceColor = ccolor;      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

plot( log10( rPIE_alpha_T( data_slice ) ), epoch_sigopt_mean( data_slice ), '-', 'Marker', '.', 'linewidth', 2, 'markersize', 30, 'color', ccolor )

hold off

grid on
set( gca, 'fontweight', 'bold' )
ylabel('Epoch')
xlabel('log10( rPIE_alpha_T )','Interpreter','none')
title( num2str( log10_sigma_optim, 'log10_sigma_optim  = %.3f'),'Interpreter','none')
% title( num2str( [ log10_sigma_optim, rPIE_alpha_phi( data_slice( end )) ], 'log10_sigma_optim  = %.3f, rPIE_alpha_phi = %.1e'),'Interpreter','none')
ylim([ 1, 6000 ])


legend( '',                      ...
        'rPIE_alpha_phi = 1e-0', ...
        '',                      ...
        'rPIE_alpha_phi = 1e-3', ...
        '',                      ...
        'rPIE_alpha_phi = 1e-6', ...
        '',                      ...
        'rPIE_alpha_phi = 1e-9', ...
        'Interpreter','none' );
    
    
% legend( 'rPIE_alpha_phi = 1e-0',         ...
%         'rPIE_alpha_phi = 1e-0 min/max', ...
%         'rPIE_alpha_phi = 1e-3',         ...
%         'rPIE_alpha_phi = 1e-3 min/max', ...
%         'rPIE_alpha_phi = 1e-6',         ...
%         'rPIE_alpha_phi = 1e-6 min/max', ...
%         'rPIE_alpha_phi = 1e-9',         ...
%         'rPIE_alpha_phi = 1e-9 min/max', ...
%         'Interpreter','none' );
%     

    
% legend( 'rPIE_alpha_phi = 1e-0',         ...
%         '', ...
%         'rPIE_alpha_phi = 1e-3',         ...
%         '', ...
%         'rPIE_alpha_phi = 1e-6',         ...
%         '', ...
%         'rPIE_alpha_phi = 1e-9',         ...
%         '', ...
%         'Interpreter','none' );

5;


%================================================


%========

h1 = figure;  
set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )

%========

data_slice = [ 1, 5, 9, 13 ];
ccolor = [ 0.0, 0.0, 0.8 ];

x     = log10( rPIE_alpha_phi( data_slice ));
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ max_epoch_sigopt_10trials( data_slice ), fliplr( min_epoch_sigopt_10trials( data_slice ) ) ];

hold on

p = fill( xconf, yconf, ccolor );
p.FaceColor = ccolor;      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

plot( log10( rPIE_alpha_phi( data_slice ) ), epoch_sigopt_mean( data_slice ), '-', 'Marker', '.', 'linewidth', 2, 'markersize', 30, 'color', ccolor )

hold off

grid on
set( gca, 'fontweight', 'bold' )
ylabel('Epoch')
xlabel('log10( rPIE_alpha_T )','Interpreter','none')

% title( num2str( [ log10_sigma_optim, rPIE_alpha_phi( data_slice( end )) ], 'log10_sigma_optim  = %.3f, rPIE_alpha_phi = %.1e'),'Interpreter','none')
title( num2str( log10_sigma_optim, 'log10_sigma_optim  = %.3f'),'Interpreter','none')
% ylim([ 1, 5000 ])

%========

data_slice = [ 1, 5, 9, 13 ] + 1;
ccolor = [ 0.0, 0.8, 0.0 ];

x     = log10( rPIE_alpha_phi( data_slice ));
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ max_epoch_sigopt_10trials( data_slice ), fliplr( min_epoch_sigopt_10trials( data_slice ) ) ];

hold on

p = fill( xconf, yconf, ccolor );
p.FaceColor = ccolor;      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

plot( log10( rPIE_alpha_phi( data_slice ) ), epoch_sigopt_mean( data_slice ), '-', 'Marker', '.', 'linewidth', 2, 'markersize', 30, 'color', ccolor )

hold off

grid on
set( gca, 'fontweight', 'bold' )
ylabel('Epoch')
xlabel('log10( rPIE_alpha_T )','Interpreter','none')
title( num2str( log10_sigma_optim, 'log10_sigma_optim  = %.3f'),'Interpreter','none')
% title( num2str( [ log10_sigma_optim, rPIE_alpha_phi( data_slice( end )) ], 'log10_sigma_optim  = %.3f, rPIE_alpha_phi = %.1e'),'Interpreter','none')
% ylim([ 1, 5000 ])

%========

data_slice = [ 1, 5, 9, 13 ] + 2;
ccolor = [ 0.8, 0.0, 0.0 ];

x     = log10( rPIE_alpha_phi( data_slice ));
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ max_epoch_sigopt_10trials( data_slice ), fliplr( min_epoch_sigopt_10trials( data_slice ) ) ];

hold on

p = fill( xconf, yconf, ccolor );
p.FaceColor = ccolor;      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

plot( log10( rPIE_alpha_phi( data_slice ) ), epoch_sigopt_mean( data_slice ), '-', 'Marker', '.', 'linewidth', 2, 'markersize', 30, 'color', ccolor )

hold off

grid on
set( gca, 'fontweight', 'bold' )
ylabel('Epoch')
xlabel('log10( rPIE_alpha_T )','Interpreter','none')
title( num2str( log10_sigma_optim, 'log10_sigma_optim  = %.3f'),'Interpreter','none')
% title( num2str( [ log10_sigma_optim, rPIE_alpha_phi( data_slice( end )) ], 'log10_sigma_optim  = %.3f, rPIE_alpha_phi = %.1e'),'Interpreter','none')
% ylim([ 1, 6000 ])

%========

data_slice = [ 1, 5, 9, 13 ] + 3;
ccolor = [ 0.8, 0.4, 0.8 ];

x     = log10( rPIE_alpha_phi( data_slice ));
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ max_epoch_sigopt_10trials( data_slice ), fliplr( min_epoch_sigopt_10trials( data_slice ) ) ];

hold on

p = fill( xconf, yconf, ccolor );
p.FaceColor = ccolor;      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;

plot( log10( rPIE_alpha_phi( data_slice ) ), epoch_sigopt_mean( data_slice ), '-', 'Marker', '.', 'linewidth', 2, 'markersize', 30, 'color', ccolor )

hold off

grid on
set( gca, 'fontweight', 'bold' )
ylabel('Epoch')
xlabel('log10( rPIE_alpha_T )','Interpreter','none')
title( num2str( log10_sigma_optim, 'log10_sigma_optim  = %.3f'),'Interpreter','none')
% title( num2str( [ log10_sigma_optim, rPIE_alpha_phi( data_slice( end )) ], 'log10_sigma_optim  = %.3f, rPIE_alpha_phi = %.1e'),'Interpreter','none')
ylim([ 1, 6000 ])



legend( '',                      ...
        'rPIE_alpha_phi = 1e-0', ...
        '',                      ...
        'rPIE_alpha_phi = 1e-3', ...
        '',                      ...
        'rPIE_alpha_phi = 1e-6', ...
        '',                      ...
        'rPIE_alpha_phi = 1e-9', ...
        'Interpreter','none' );
    
    





































return





%========

kk = 1;

for ii = 1 : 4
    
    epoch_indx = find( metrics.mb20_noise_rmbg{ ii }.log10_meas_all_avg < log10_sigma_optim, 1 );
    epoch_sigopt_mean( kk ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( epoch_indx );
    rPIE_alpha_T( kk )      = metrics.mb20_noise_rmbg{ ii }.rPIE_alpha_T( 1 );   
    rPIE_alpha_phi( kk )    = metrics.mb20_noise_rmbg{ ii }.rPIE_alpha_phi( 1 );    

    for cc = 1 : 10

        II = find( metrics.mb20_noise_rmbg{ ii }.log10_meas_all( :, cc ) < log10_sigma_optim, 1 );

        if isempty( II )

            tmp0( cc ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( end );

        else

            tmp0( cc ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( II );

        end

    end

    min_epoch_sigopt_10trials( kk ) = min( tmp0 );
    max_epoch_sigopt_10trials( kk ) = max( tmp0 );
    
    kk = kk + 1;
    
end

ccolor = [ 0.0, 0.0, 0.8 ];
% ccolor = [ 0.0, 0.8, 0.0 ];
% ccolor = [ 0.8, 0.0, 0.0 ];

h1 = figure;  
set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )

x     = log10( rPIE_alpha_T );
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ max_epoch_sigopt_10trials, fliplr( min_epoch_sigopt_10trials ) ];

hold on
p = fill( xconf, yconf, ccolor );
p.FaceColor = ccolor;      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.1;
hold off

hold on
plot( log10( rPIE_alpha_T ), epoch_sigopt_mean, '-', 'Marker', '.', 'linewidth', 2, 'markersize', 30, 'color', ccolor )
hold off

grid on
set( gca, 'fontweight', 'bold' )
ylabel('Epoch')
xlabel('log10( rPIE_alpha_T )','Interpreter','none')

title( num2str( [ log10_sigma_optim, rPIE_alpha_phi( 1 ) ], 'log10_sigma_optim  = %.3f, rPIE_alpha_phi = %.1e'),'Interpreter','none')
ylim([ 1, 5000 ])

clear( 'rPIE_alpha_phi', 'epoch_sigopt_mean', 'rPIE_alpha_T', 'max_epoch_sigopt_10trials', 'min_epoch_sigopt_10trials' )

%========

kk = 1;

for ii = 5 : 8

    epoch_indx = find( metrics.mb20_noise_rmbg{ ii }.log10_meas_all_avg < log10_sigma_optim, 1 );
    epoch_sigopt_mean( kk ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( epoch_indx );
    rPIE_alpha_T( kk )      = metrics.mb20_noise_rmbg{ ii }.rPIE_alpha_T( 1 );   
    rPIE_alpha_phi( kk )    = metrics.mb20_noise_rmbg{ ii }.rPIE_alpha_phi( 1 );    

    for cc = 1 : 10

        II = find( metrics.mb20_noise_rmbg{ ii }.log10_meas_all( :, cc ) < log10_sigma_optim, 1 );

        if isempty( II )

            tmp0( cc ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( end );

        else

            tmp0( cc ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( II );

        end

    end

    min_epoch_sigopt_10trials( kk ) = min( tmp0 );
    max_epoch_sigopt_10trials( kk ) = max( tmp0 );
    
    kk = kk + 1;
    
end

ccolor = [ 0.0, 0.8, 0.0 ];

h1 = figure;  
set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )

x     = log10( rPIE_alpha_T );
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ max_epoch_sigopt_10trials, fliplr( min_epoch_sigopt_10trials ) ];

hold on
p = fill( xconf, yconf, ccolor );
p.FaceColor = ccolor;      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.1;
hold off

hold on
plot( log10( rPIE_alpha_T ), epoch_sigopt_mean, '-', 'Marker', '.', 'linewidth', 2, 'markersize', 30, 'color', ccolor )
hold off

grid on
set( gca, 'fontweight', 'bold' )
ylabel('Epoch')
xlabel('log10( rPIE_alpha_T )','Interpreter','none')

title( num2str( [ log10_sigma_optim, rPIE_alpha_phi( 1 ) ], 'log10_sigma_optim  = %.3f, rPIE_alpha_phi = %.1e'),'Interpreter','none')
ylim([ 1, 5000 ])

clear( 'rPIE_alpha_phi', 'epoch_sigopt_mean', 'rPIE_alpha_T', 'max_epoch_sigopt_10trials', 'min_epoch_sigopt_10trials' )

%========

kk = 1;

for ii = 9 : 12
 
    epoch_indx = find( metrics.mb20_noise_rmbg{ ii }.log10_meas_all_avg < log10_sigma_optim, 1 );
    epoch_sigopt_mean( kk ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( epoch_indx );
    rPIE_alpha_T( kk )      = metrics.mb20_noise_rmbg{ ii }.rPIE_alpha_T( 1 );   
    rPIE_alpha_phi( kk )    = metrics.mb20_noise_rmbg{ ii }.rPIE_alpha_phi( 1 );    

    for cc = 1 : 10

        II = find( metrics.mb20_noise_rmbg{ ii }.log10_meas_all( :, cc ) < log10_sigma_optim, 1 );

        if isempty( II )

            tmp0( cc ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( end );

        else

            tmp0( cc ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( II );

        end

    end

    min_epoch_sigopt_10trials( kk ) = min( tmp0 );
    max_epoch_sigopt_10trials( kk ) = max( tmp0 );
    
    kk = kk + 1;
    
end

ccolor = [ 0.8, 0.0, 0.0 ];

h1 = figure;  
set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )

x     = log10( rPIE_alpha_T );
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ max_epoch_sigopt_10trials, fliplr( min_epoch_sigopt_10trials ) ];

hold on
p = fill( xconf, yconf, ccolor );
p.FaceColor = ccolor;      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;
hold off

hold on
plot( log10( rPIE_alpha_T ), epoch_sigopt_mean, '-', 'Marker', '.', 'linewidth', 2, 'markersize', 30, 'color', ccolor )
hold off

grid on
set( gca, 'fontweight', 'bold' )
ylabel('Epoch')
xlabel('log10( rPIE_alpha_T )','Interpreter','none')

title( num2str( [ log10_sigma_optim, rPIE_alpha_phi( 1 ) ], 'log10_sigma_optim  = %.3f, rPIE_alpha_phi = %.1e'),'Interpreter','none')
ylim([ 1, 5000 ])

clear( 'rPIE_alpha_phi', 'epoch_sigopt_mean', 'rPIE_alpha_T', 'max_epoch_sigopt_10trials', 'min_epoch_sigopt_10trials' )

%========

kk = 1;

for ii = 13 : 16
    
    epoch_indx = find( metrics.mb20_noise_rmbg{ ii }.log10_meas_all_avg < log10_sigma_optim, 1 );
    epoch_sigopt_mean( kk ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( epoch_indx );
    rPIE_alpha_T( kk )      = metrics.mb20_noise_rmbg{ ii }.rPIE_alpha_T( 1 );   
    rPIE_alpha_phi( kk )    = metrics.mb20_noise_rmbg{ ii }.rPIE_alpha_phi( 1 );    

    for cc = 1 : 10

        II = find( metrics.mb20_noise_rmbg{ ii }.log10_meas_all( :, cc ) < log10_sigma_optim, 1 );

        if isempty( II )

            tmp0( cc ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( end );

        else

            tmp0( cc ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( II );

        end

    end

    min_epoch_sigopt_10trials( kk ) = min( tmp0 );
    max_epoch_sigopt_10trials( kk ) = max( tmp0 );
    
    kk = kk + 1;
    
end

ccolor = [ 0.0, 0.0, 0.0 ];

h1 = figure;  
set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )

x     = log10( rPIE_alpha_T );
xconf = [ x, x( end : -1 : 1) ]; 
yconf = [ max_epoch_sigopt_10trials, fliplr( min_epoch_sigopt_10trials ) ];

hold on
p = fill( xconf, yconf, ccolor );
p.FaceColor = ccolor;      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.1;
hold off

hold on
plot( log10( rPIE_alpha_T ), epoch_sigopt_mean, '-', 'Marker', '.', 'linewidth', 2, 'markersize', 30, 'color', ccolor )
hold off

grid on
set( gca, 'fontweight', 'bold' )
ylabel('Epoch')
xlabel('log10( rPIE_alpha_T )','Interpreter','none')

title( num2str( [ log10_sigma_optim, rPIE_alpha_phi( 1 ) ], 'log10_sigma_optim  = %.3f, rPIE_alpha_phi = %.1e'),'Interpreter','none')
ylim([ 1, 5000 ])

clear( 'rPIE_alpha_phi', 'epoch_sigopt_mean', 'rPIE_alpha_T', 'max_epoch_sigopt_10trials', 'min_epoch_sigopt_10trials' )

end





% 
% figure
% 
% x1 = 0:0.1:40;
% y1 = 4.*cos(x1)./(x1+2);
% x2 = 1:0.2:20;
% y2 = x2.^2./x2.^3;
% 
% t = tiledlayout(1,1);
% 
% ax1 = axes(t);
% plot(ax1,x1,y1,'-r')
% 
% % ax1.Color = 'none';
% ax1.XColor = [ 0.8, 0.0, 0.0 ];
% ax1.YColor = [ 0.8, 0.0, 0.0 ];
% 
% ax2 = axes(t);
% plot(ax2,x2,y2,'-k')
% 
% ax2.XAxisLocation = 'top';
% ax2.YAxisLocation = 'right';
% ax2.YLabel.String = 'AAA';
% ax2.XLabel.String = 'BBB';
% 
% 
% ax2.Color = 'none';
% % ax2.Color = [ 0.0, 0.7, 0.0 ];
% ax2.XColor = [ 0.0, 0.7, 0.0 ];
% ax2.YColor = [ 0.0, 0.7, 0.0 ];
% 
% ax2.XLabel.Color = [ 0.0, 0.0, 0.5 ];
% 
% ax1.Box = 'off';
% ax2.Box = 'off';
% 
% 
% 
% 
% figure
% x1 = 0:0.1:40;
% y1 = 4.*cos(x1)./(x1+2);
% x2 = fliplr( x1 );
% 
% t = tiledlayout(1,1);
% 
% ax1 = axes(t);
% plot(ax1,x1,y1,'-r')
% set( ax1, 'xlabel', 'zzz' )
% 
% 
% ax2 = axes(t);
% plot(ax2,x2,y1,'-k')
% ax2.XAxisLocation = 'top';
% set( ax2, 'xlabel', 'hhh' )


%{

figure

x = 1:4;
y = rand(1,4);
plot(x,y);

ax = gca;

% ax.XTick = [1 2 3 4];
ax.XTick = x;

ax.XTickLabel = '';

myLabels = { '1', '2', '3', '4'; 
             'Line2a', 'Line2b', 'Line2c', 'Line2d';
             'Line3a', 'Line3b', 'Line3c', 'Line3d'; 
             'Line4a', 'Line4b', 'Line4c', 'Line4d'; 
             'Line5a', 'Line5b', 'Line5c', 'Line5d' };

% for i = 1:length(myLabels)
%     text(i, ax.YLim(1), sprintf('%s\n%s\n%s', myLabels{:,i}), ...
%         'horizontalalignment', 'center', 'verticalalignment', 'top');    
% end

ccolors = [ [ 0.8, 0.0, 0.0 ]; [ 0.0, 0.8, 0.0 ]; [ 0.0, 0.0, 0.8 ]; [ 0.0, 0.0, 0.0 ]; [ 0.8, 0.8, 0.0 ] ];

for i = 1 : size( myLabels, 2 )
    
    text(ax.XTick(i), ax.YLim(1), sprintf('%s\n%s\n%s\n%s\n%s', myLabels{:,i}), 'horizontalalignment', 'center', 'verticalalignment', 'top' );    
    
end

ax.XLabel.String = sprintf('\n\n\n\n\n%s', 'X-Axis Label');



5;


























return



figure
% yyaxis left
x = linspace(0,10);
y = sin(3*x);
plot(x,y)
yyaxis right
y2 = sin(3*x).*exp(0.5*x);
plot(x,y2)



figure

x1 = 0:0.1:40;
y1 = 4.*cos(x1)./(x1+2);
x2 = 1:0.2:20;
y2 = x2.^2./x2.^3;

t = tiledlayout(1,1);

ax1 = axes(t);
plot(ax1,x1,y1,'-r')

% ax1.Color = 'none';
ax1.XColor = [ 0.8, 0.0, 0.0 ];
ax1.YColor = [ 0.8, 0.0, 0.0 ];

ax2 = axes(t);
plot(ax2,x2,y2,'-k')

ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';

ax2.Color = 'none';
% ax2.Color = [ 0.0, 0.7, 0.0 ];
ax2.XColor = [ 0.0, 0.7, 0.0 ];
ax2.YColor = [ 0.0, 0.7, 0.0 ];

ax1.Box = 'off';
ax2.Box = 'off';














%}








%====================================================================================================================================================
% 
% function duchii_plots_variance( metrics )
% 
% 
% log10_sigma_optim = 4.1;
% 
% %========
% 
% kk = 1;
% 
% for ii = 1 : 4
%     
%     epoch_indx = find( metrics.mb20_noise_rmbg{ ii }.log10_meas_all_avg < log10_sigma_optim, 1 );
%     epoch_sigopt_mean( kk ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( epoch_indx );
%     rPIE_alpha_T( kk )      = metrics.mb20_noise_rmbg{ ii }.rPIE_alpha_T( 1 );   
%     rPIE_alpha_phi( kk )    = metrics.mb20_noise_rmbg{ ii }.rPIE_alpha_phi( 1 );    
% 
%     for cc = 1 : 10
% 
%         II = find( metrics.mb20_noise_rmbg{ ii }.log10_meas_all( :, cc ) < log10_sigma_optim, 1 );
% 
%         if isempty( II )
% 
%             tmp0( cc ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( end );
% 
%         else
% 
%             tmp0( cc ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( II );
% 
%         end
% 
%     end
% 
%     variance_epoch_sigopt_10trials( kk ) = std( tmp0 );
%     
%     kk = kk + 1;
%     
% end
% 
% ccolor = [ 0.0, 0.0, 0.8 ];
% % ccolor = [ 0.0, 0.8, 0.0 ];
% % ccolor = [ 0.8, 0.0, 0.0 ];
% 
% h1 = figure;  
% set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
% 
% x     = log10( rPIE_alpha_T );
% xconf = [ x, x( end : -1 : 1) ]; 
% yconf = [ epoch_sigopt_mean + 0.5 * variance_epoch_sigopt_10trials, fliplr( epoch_sigopt_mean - 0.5 * variance_epoch_sigopt_10trials ) ];
% 
% hold on
% p = fill( xconf, yconf, ccolor );
% p.FaceColor = ccolor;      
% p.EdgeColor = 'none';           
% p.FaceAlpha = 0.1;
% hold off
% 
% hold on
% plot( log10( rPIE_alpha_T ), epoch_sigopt_mean, '-', 'Marker', '.', 'linewidth', 2, 'markersize', 30, 'color', ccolor )
% hold off
% 
% grid on
% set( gca, 'fontweight', 'bold' )
% ylabel('Epoch')
% xlabel('log10( rPIE_alpha_T )','Interpreter','none')
% 
% title( num2str( [ log10_sigma_optim, rPIE_alpha_phi( 1 ) ], 'log10_sigma_optim  = %.3f, rPIE_alpha_phi = %.1e'),'Interpreter','none')
% ylim([ 1, 5000 ])
% 
% clear( 'rPIE_alpha_phi', 'epoch_sigopt_mean', 'rPIE_alpha_T', 'max_epoch_sigopt_10trials', 'min_epoch_sigopt_10trials' )
% 
% %========
% 
% kk = 1;
% 
% for ii = 5 : 8
% 
%     epoch_indx = find( metrics.mb20_noise_rmbg{ ii }.log10_meas_all_avg < log10_sigma_optim, 1 );
%     epoch_sigopt_mean( kk ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( epoch_indx );
%     rPIE_alpha_T( kk )      = metrics.mb20_noise_rmbg{ ii }.rPIE_alpha_T( 1 );   
%     rPIE_alpha_phi( kk )    = metrics.mb20_noise_rmbg{ ii }.rPIE_alpha_phi( 1 );    
% 
%     for cc = 1 : 10
% 
%         II = find( metrics.mb20_noise_rmbg{ ii }.log10_meas_all( :, cc ) < log10_sigma_optim, 1 );
% 
%         if isempty( II )
% 
%             tmp0( cc ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( end );
% 
%         else
% 
%             tmp0( cc ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( II );
% 
%         end
% 
%     end
% 
%     min_epoch_sigopt_10trials( kk ) = min( tmp0 );
%     max_epoch_sigopt_10trials( kk ) = max( tmp0 );
%     
%     kk = kk + 1;
%     
% end
% 
% ccolor = [ 0.0, 0.8, 0.0 ];
% 
% h1 = figure;  
% set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
% 
% x     = log10( rPIE_alpha_T );
% xconf = [ x, x( end : -1 : 1) ]; 
% yconf = [ max_epoch_sigopt_10trials, fliplr( min_epoch_sigopt_10trials ) ];
% 
% hold on
% p = fill( xconf, yconf, ccolor );
% p.FaceColor = ccolor;      
% p.EdgeColor = 'none';           
% p.FaceAlpha = 0.1;
% hold off
% 
% hold on
% plot( log10( rPIE_alpha_T ), epoch_sigopt_mean, '-', 'Marker', '.', 'linewidth', 2, 'markersize', 30, 'color', ccolor )
% hold off
% 
% grid on
% set( gca, 'fontweight', 'bold' )
% ylabel('Epoch')
% xlabel('log10( rPIE_alpha_T )','Interpreter','none')
% 
% title( num2str( [ log10_sigma_optim, rPIE_alpha_phi( 1 ) ], 'log10_sigma_optim  = %.3f, rPIE_alpha_phi = %.1e'),'Interpreter','none')
% ylim([ 1, 5000 ])
% 
% clear( 'rPIE_alpha_phi', 'epoch_sigopt_mean', 'rPIE_alpha_T', 'max_epoch_sigopt_10trials', 'min_epoch_sigopt_10trials' )
% 
% %========
% 
% kk = 1;
% 
% for ii = 9 : 12
%  
%     epoch_indx = find( metrics.mb20_noise_rmbg{ ii }.log10_meas_all_avg < log10_sigma_optim, 1 );
%     epoch_sigopt_mean( kk ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( epoch_indx );
%     rPIE_alpha_T( kk )      = metrics.mb20_noise_rmbg{ ii }.rPIE_alpha_T( 1 );   
%     rPIE_alpha_phi( kk )    = metrics.mb20_noise_rmbg{ ii }.rPIE_alpha_phi( 1 );    
% 
%     for cc = 1 : 10
% 
%         II = find( metrics.mb20_noise_rmbg{ ii }.log10_meas_all( :, cc ) < log10_sigma_optim, 1 );
% 
%         if isempty( II )
% 
%             tmp0( cc ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( end );
% 
%         else
% 
%             tmp0( cc ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( II );
% 
%         end
% 
%     end
% 
%     min_epoch_sigopt_10trials( kk ) = min( tmp0 );
%     max_epoch_sigopt_10trials( kk ) = max( tmp0 );
%     
%     kk = kk + 1;
%     
% end
% 
% ccolor = [ 0.8, 0.0, 0.0 ];
% 
% h1 = figure;  
% set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
% 
% x     = log10( rPIE_alpha_T );
% xconf = [ x, x( end : -1 : 1) ]; 
% yconf = [ max_epoch_sigopt_10trials, fliplr( min_epoch_sigopt_10trials ) ];
% 
% hold on
% p = fill( xconf, yconf, ccolor );
% p.FaceColor = ccolor;      
% p.EdgeColor = 'none';           
% p.FaceAlpha = 0.2;
% hold off
% 
% hold on
% plot( log10( rPIE_alpha_T ), epoch_sigopt_mean, '-', 'Marker', '.', 'linewidth', 2, 'markersize', 30, 'color', ccolor )
% hold off
% 
% grid on
% set( gca, 'fontweight', 'bold' )
% ylabel('Epoch')
% xlabel('log10( rPIE_alpha_T )','Interpreter','none')
% 
% title( num2str( [ log10_sigma_optim, rPIE_alpha_phi( 1 ) ], 'log10_sigma_optim  = %.3f, rPIE_alpha_phi = %.1e'),'Interpreter','none')
% ylim([ 1, 5000 ])
% 
% clear( 'rPIE_alpha_phi', 'epoch_sigopt_mean', 'rPIE_alpha_T', 'max_epoch_sigopt_10trials', 'min_epoch_sigopt_10trials' )
% 
% %========
% 
% kk = 1;
% 
% for ii = 13 : 16
%     
%     epoch_indx = find( metrics.mb20_noise_rmbg{ ii }.log10_meas_all_avg < log10_sigma_optim, 1 );
%     epoch_sigopt_mean( kk ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( epoch_indx );
%     rPIE_alpha_T( kk )      = metrics.mb20_noise_rmbg{ ii }.rPIE_alpha_T( 1 );   
%     rPIE_alpha_phi( kk )    = metrics.mb20_noise_rmbg{ ii }.rPIE_alpha_phi( 1 );    
% 
%     for cc = 1 : 10
% 
%         II = find( metrics.mb20_noise_rmbg{ ii }.log10_meas_all( :, cc ) < log10_sigma_optim, 1 );
% 
%         if isempty( II )
% 
%             tmp0( cc ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( end );
% 
%         else
% 
%             tmp0( cc ) = metrics.mb20_noise_rmbg{ ii }.it( 1 ).mtot( II );
% 
%         end
% 
%     end
% 
%     min_epoch_sigopt_10trials( kk ) = min( tmp0 );
%     max_epoch_sigopt_10trials( kk ) = max( tmp0 );
%     
%     kk = kk + 1;
%     
% end
% 
% ccolor = [ 0.0, 0.0, 0.0 ];
% 
% h1 = figure;  
% set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
% 
% x     = log10( rPIE_alpha_T );
% xconf = [ x, x( end : -1 : 1) ]; 
% yconf = [ max_epoch_sigopt_10trials, fliplr( min_epoch_sigopt_10trials ) ];
% 
% hold on
% p = fill( xconf, yconf, ccolor );
% p.FaceColor = ccolor;      
% p.EdgeColor = 'none';           
% p.FaceAlpha = 0.1;
% hold off
% 
% hold on
% plot( log10( rPIE_alpha_T ), epoch_sigopt_mean, '-', 'Marker', '.', 'linewidth', 2, 'markersize', 30, 'color', ccolor )
% hold off
% 
% grid on
% set( gca, 'fontweight', 'bold' )
% ylabel('Epoch')
% xlabel('log10( rPIE_alpha_T )','Interpreter','none')
% 
% title( num2str( [ log10_sigma_optim, rPIE_alpha_phi( 1 ) ], 'log10_sigma_optim  = %.3f, rPIE_alpha_phi = %.1e'),'Interpreter','none')
% ylim([ 1, 5000 ])
% 
% clear( 'rPIE_alpha_phi', 'epoch_sigopt_mean', 'rPIE_alpha_T', 'max_epoch_sigopt_10trials', 'min_epoch_sigopt_10trials' )
% 
% 
% 
% end