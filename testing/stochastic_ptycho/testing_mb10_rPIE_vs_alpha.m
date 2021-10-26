% 

%
%{

clear; close all; testing_mb10_rPIE_vs_alpha

%}

%====================================================================================================================================================

addpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/shadedErrorBar/' );  

%====================================================================================================================================================

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/mb0p10/';

%====================================================================================================================================================

path_data = {};
N_trials  = [];

%========

% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha1p00/independenttrials_16Aug2021_t024241/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p80/independenttrials_15Aug2021_t115003/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p50/independenttrials_14Aug2021_t210351/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p30/independenttrials_14Aug2021_t061851/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p10/independenttrials_12Aug2021_t125535/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p09/independenttrials_12Aug2021_t044127/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p08/independenttrials_11Aug2021_t202809/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p07/independenttrials_11Aug2021_t121230/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p06/independenttrials_11Aug2021_t035951/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p05/independenttrials_10Aug2021_t195103/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p04/independenttrials_10Aug2021_t113245/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p03/independenttrials_10Aug2021_t031506/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p02/independenttrials_09Aug2021_t164930/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p01_a/independenttrials_13Aug2021_t154650/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p01_b/independenttrials_06Aug2021_t095015/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p0001/independenttrials_28Aug2021_t145156/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p00001/independenttrials_27Aug2021_t131958/' ];
% N_trials( end + 1 )  = 10;

%========

##

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha1p00_oldsample_for_probe_update/independenttrials_23Oct2021_t190740/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p80/independenttrials_15Aug2021_t115003/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p50/independenttrials_14Aug2021_t210351/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p30/independenttrials_14Aug2021_t061851/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p10/independenttrials_12Aug2021_t125535/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p09/independenttrials_12Aug2021_t044127/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p08/independenttrials_11Aug2021_t202809/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p07/independenttrials_11Aug2021_t121230/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p06/independenttrials_11Aug2021_t035951/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p05/independenttrials_10Aug2021_t195103/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p04/independenttrials_10Aug2021_t113245/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p03/independenttrials_10Aug2021_t031506/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p02/independenttrials_09Aug2021_t164930/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p01_a/independenttrials_13Aug2021_t154650/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p01_b/independenttrials_06Aug2021_t095015/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p0001/independenttrials_28Aug2021_t145156/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p00001/independenttrials_27Aug2021_t131958/' ];
N_trials( end + 1 )  = 10;

%====================================================================================================================================================

for jj = 1 : length( path_data )
    
    metrics( jj ) = load_and_plot( path_data{ jj }, N_trials( jj ) ); %#ok<SAGROW>

end

%====================================================================================================================================================

function metrics = load_and_plot( path_data, N_trials )

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

    name_data = num2str( [ sim_ptycho2DTPA{ 1 }.sol.rPIE_alpha, sim_ptycho2DTPA{ 1 }.sol.spos.rand_spos_subset_pct ], 'rPIE_alpha = %0.8f, MBpct = %0.4f');


    %=================================================
    % create a struct from metrics data for future use
    %=================================================
    
    for ii = 1 : N_trials

        metrics.meas_all( :, ii )          = sim_ptycho2DTPA{ ii }.sol.metrics.meas_all;
        metrics.timing( ii )               = sim_ptycho2DTPA{ ii }.sol.timings;
        metrics.it( ii )                   = sim_ptycho2DTPA{ ii }.sol.it;
        metrics.rPIE_alpha( ii )           = sim_ptycho2DTPA{ ii }.sol.rPIE_alpha;
        metrics.rand_spos_subset_pct( ii ) = sim_ptycho2DTPA{ ii }.sol.spos.rand_spos_subset_pct;
        
    end
    
    name_data = num2str( [ metrics.rPIE_alpha( 1 ), metrics.rand_spos_subset_pct( 1 ) ], 'rPIE_alpha = %0.8f, MBpct = %0.4f');
    
    metrics.meas_all_avg_of_log = meas_all_avg;
    metrics.name_data           = name_data;

    %========
    
    skip = 1;
    
    h1 = figure();  
    set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )

    x     = metrics.it( 1 ).mtot( 1 : skip : end );
    xconf = [ x, x( end : -1 : 1) ]; 
    yconf = [ log10_meas_all_max, fliplr( log10_meas_all_min ) ];
    
    p = fill( xconf, yconf, 'red' );
    p.FaceColor = [ 1, 0.8, 0.8 ];      
    p.EdgeColor = 'none';           

    hold on

    plot( metrics.it( 1 ).mtot( 1 : skip : end ), meas_all_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )

     
    xlabel('Epoch')
    ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
    hold off
    title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
    grid on
    ylim( y_lim )
%     legend( 'Location', 'northeast' ) 
    
    hold off
     
    %========
    
%     skip = 1;
% 
%     h1 = figure();  
%     set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
% 
%     hold on
% 
%     for ii = 1 : N_trials
% 
%         plot( metrics.it( ii ).mtot( 1 : skip : end ), log10( metrics.meas_all( 1 : skip : end, ii )  ), '-', 'linewidth', 2 )
% 
%     end
% 
%     plot( metrics.it( 1 ).mtot( 1 : skip : end ), meas_all_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
%     xlabel('Epoch')
%     ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
%     hold off
%     title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
%     grid on
%     ylim( y_lim )
%     legend( 'Location', 'northeast' ) 
%     
%     hold off

end





























% %========
% 
% for jj = 1 : length( path_data )
%     
%     load_and_plot( path_data{ jj }, N_trials( jj ) )
%     
% end
% 
% %====================================================================================================================================================
% 
% function load_and_plot( path_data, N_trials )
% 
%     y_lim = [-1, 6];
%     
%     %========
%     
%     sim_ptycho2DTPA = cell( N_trials, 1 );
%     
%     for ii = 1 : N_trials
% 
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
%     end
% 
%     %========
% 
%     meas_all_avg = 0;
% 
%     for ii = 1 : N_trials
% 
%         meas_all_avg = meas_all_avg + log10( sim_ptycho2DTPA{ ii }.sol.metrics.meas_all ); 
% 
%     end
% 
%     meas_all_avg = meas_all_avg / N_trials;
% 
%     %========
% 
%     skip = 1;
% 
% %     figure; 
%     h1 = figure();  
%     set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
% 
%     hold on
% 
%     for ii = 1 : N_trials
% 
%         plot( sim_ptycho2DTPA{ ii }.sol.it.mtot( 1 : skip : end ), log10( sim_ptycho2DTPA{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
%     end
% 
%     plot( sim_ptycho2DTPA{ ii }.sol.it.mtot( 1 : skip : end ), meas_all_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
%     name_data = num2str( [ sim_ptycho2DTPA{ii}.sol.rPIE_alpha, sim_ptycho2DTPA{ii}.sol.spos.rand_spos_subset_pct ], 'rPIE_alpha = %0.8f, MBpct = %0.4f');
%     
%     xlabel('Epoch')
%     ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
%     hold off
%     title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
%     grid on
%     ylim( y_lim )
%     legend( 'Location', 'northeast' ) 
% 
% 
% end
% 
% %====================================================================================================================================================

















































% 
% 
% 
% 
% 
% y_lim = [-1, 6];
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha01_gpu2/independenttrials_13Aug2021_t154650/';
%     cdi_ER_rPIE_0p01{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p01_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p01_avg = cdi_ER_rPIE_0p01_avg + log10( cdi_ER_rPIE_0p01{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p01_avg = cdi_ER_rPIE_0p01_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p01{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p01{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p01_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.01, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear( 'cdi_ER_rPIE_0p01' )
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha01_gpu3/independenttrials_06Aug2021_t095015/';
%     cdi_ER_rPIE_0p01{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p01_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p01_avg = cdi_ER_rPIE_0p01_avg + log10( cdi_ER_rPIE_0p01{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p01_avg = cdi_ER_rPIE_0p01_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p01{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p01{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p01_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.01, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear( 'cdi_ER_rPIE_0p01' )
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha02_gpu2/independenttrials_09Aug2021_t164930/';
%     cdi_ER_rPIE_0p02{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p02_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p02_avg = cdi_ER_rPIE_0p02_avg + log10( cdi_ER_rPIE_0p02{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p02_avg = cdi_ER_rPIE_0p02_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p02{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p02{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p02{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p02_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.02, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear( 'cdi_ER_rPIE_0p02' )
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha03_gpu2/independenttrials_10Aug2021_t031506/';
%     cdi_ER_rPIE_0p03{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p03_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p03_avg = cdi_ER_rPIE_0p03_avg + log10( cdi_ER_rPIE_0p03{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p03_avg = cdi_ER_rPIE_0p03_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p03{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p03{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p03{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p03_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.03, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear( 'cdi_ER_rPIE_0p03' )
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha04_gpu2/independenttrials_10Aug2021_t113245/';
%     cdi_ER_rPIE_0p04{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p04_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p04_avg = cdi_ER_rPIE_0p04_avg + log10( cdi_ER_rPIE_0p04{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p04_avg = cdi_ER_rPIE_0p04_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p04{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p04{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p04{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p04_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.04, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear( 'cdi_ER_rPIE_0p04' )
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha05_gpu2/independenttrials_10Aug2021_t195103/';
%     cdi_ER_rPIE_0p05{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p05_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p05_avg = cdi_ER_rPIE_0p05_avg + log10( cdi_ER_rPIE_0p05{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p05_avg = cdi_ER_rPIE_0p05_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p05{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p05{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p05{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p05_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.05, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear( 'cdi_ER_rPIE_0p05' )
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha06_gpu2/independenttrials_11Aug2021_t035951/';
%     cdi_ER_rPIE_0p06{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p06_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p06_avg = cdi_ER_rPIE_0p06_avg + log10( cdi_ER_rPIE_0p06{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p06_avg = cdi_ER_rPIE_0p06_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p06{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p06{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p06{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p06_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.06, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear( 'cdi_ER_rPIE_0p06' )
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha07_gpu2/independenttrials_11Aug2021_t121230/';
%     cdi_ER_rPIE_0p07{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p07_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p07_avg = cdi_ER_rPIE_0p07_avg + log10( cdi_ER_rPIE_0p07{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p07_avg = cdi_ER_rPIE_0p07_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p07{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p07{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p07{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p07_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.07, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear( 'cdi_ER_rPIE_0p07' )
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha08_gpu2/independenttrials_11Aug2021_t202809/';
%     cdi_ER_rPIE_0p08{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p08_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p08_avg = cdi_ER_rPIE_0p08_avg + log10( cdi_ER_rPIE_0p08{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p08_avg = cdi_ER_rPIE_0p08_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p08{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p08{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p08{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p08_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.08, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear( 'cdi_ER_rPIE_0p08' )
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha09_gpu2/independenttrials_12Aug2021_t044127/';
%     cdi_ER_rPIE_0p09{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p09_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p09_avg = cdi_ER_rPIE_0p09_avg + log10( cdi_ER_rPIE_0p09{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p09_avg = cdi_ER_rPIE_0p09_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p09{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p09{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p09{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p09_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.09, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear( 'cdi_ER_rPIE_0p09' )
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha10_gpu2/independenttrials_12Aug2021_t125535/';
%     cdi_ER_rPIE_0p10{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p10_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p10_avg = cdi_ER_rPIE_0p10_avg + log10( cdi_ER_rPIE_0p10{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p10_avg = cdi_ER_rPIE_0p10_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p10{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p10{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p10{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p10_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.10, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear( 'cdi_ER_rPIE_0p10' )
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha30_gpu2/independenttrials_14Aug2021_t061851/';
%     cdi_ER_rPIE_0p30{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p30_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p30_avg = cdi_ER_rPIE_0p30_avg + log10( cdi_ER_rPIE_0p30{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p30_avg = cdi_ER_rPIE_0p30_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p30{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p30{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p30{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p30_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.30, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear( 'cdi_ER_rPIE_0p30' )
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha50_gpu2/independenttrials_14Aug2021_t210351/';
%     cdi_ER_rPIE_0p50{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p50_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p50_avg = cdi_ER_rPIE_0p50_avg + log10( cdi_ER_rPIE_0p50{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p50_avg = cdi_ER_rPIE_0p50_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p50{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p50{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p50{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p50_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.50, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear( 'cdi_ER_rPIE_0p50' )
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha80_gpu2/independenttrials_15Aug2021_t115003/';
%     cdi_ER_rPIE_0p80{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p80_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p80_avg = cdi_ER_rPIE_0p80_avg + log10( cdi_ER_rPIE_0p80{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p80_avg = cdi_ER_rPIE_0p80_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p80{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p80{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p80{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p80_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.80, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear( 'cdi_ER_rPIE_0p80' )
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha100/independenttrials_16Aug2021_t024241/';
%     cdi_ER_rPIE_0p100{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p100_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p100_avg = cdi_ER_rPIE_0p100_avg + log10( cdi_ER_rPIE_0p100{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p100_avg = cdi_ER_rPIE_0p100_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p100{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p100{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p100{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p100_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 1.00, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear( 'cdi_ER_rPIE_0p100' )
% 
% %====================================================================================================================================================
% 
% 
% 
% 
% 
% 
