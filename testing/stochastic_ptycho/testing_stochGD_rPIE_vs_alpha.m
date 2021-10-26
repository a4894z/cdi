% 

%
%{

close all; clear; testing_stochGD_rPIE_vs_alpha

%}

%====================================================================================================================================================

addpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/shadedErrorBar/' );  

%====================================================================================================================================================

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/blind_block_stoch/';

%====================================================================================================================================================

path_data = {};
N_trials  = [];

%========

% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p90/independenttrials_23Aug2021_t065834/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p80/independenttrials_22Aug2021_t125306/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p70/independenttrials_21Aug2021_t151804/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p60/independenttrials_20Aug2021_t170854/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p50/independenttrials_23Aug2021_t060643/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p40/independenttrials_22Aug2021_t121326/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p30/independenttrials_21Aug2021_t150513/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p20/independenttrials_20Aug2021_t171816/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p10/independenttrials_12Aug2021_t111313/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p09/independenttrials_11Aug2021_t194747/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p08/independenttrials_11Aug2021_t042043/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p07/independenttrials_10Aug2021_t130027/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p06/independenttrials_09Aug2021_t214746/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p05/independenttrials_06Aug2021_t093611/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p04/independenttrials_05Aug2021_t105850/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p03/independenttrials_05Aug2021_t100113/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p02/independenttrials_05Aug2021_t095936/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p01/independenttrials_04Aug2021_t150026/' ];
% N_trials( end + 1 )  = 10;

% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p005/independenttrials_14Sep2021_t113015/' ];
% N_trials( end + 1 )  = 10;

% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p001/independenttrials_13Sep2021_t164005/' ];
% N_trials( end + 1 )  = 10;

%========

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha1p00_oldsample_for_probe_update/independenttrials_10Oct2021_t140518/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p80_oldsample_for_probe_update/independenttrials_09Oct2021_t045024/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p60_oldsample_for_probe_update/independenttrials_07Oct2021_t193535/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p40_oldsample_for_probe_update/independenttrials_06Oct2021_t102502/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p20_oldsample_for_probe_update/independenttrials_05Oct2021_t011944/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p10_oldsample_for_probe_update/independenttrials_03Oct2021_t161143/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p05_oldsample_for_probe_update/independenttrials_16Sep2021_t151721/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p01_oldsample_for_probe_update/independenttrials_22Sep2021_t104744/' ];
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

    %=================================================
    % create a struct from metrics data for future use
    %=================================================
    
    for ii = 1 : N_trials

        metrics.meas_all( :, ii )          = sim_ptycho2DTPA{ ii }.sol.metrics.meas_all;
        metrics.timing( ii )               = sim_ptycho2DTPA{ ii }.sol.timings;
        metrics.it( ii )                   = sim_ptycho2DTPA{ ii }.sol.it;
        metrics.rPIE_alpha( ii )           = sim_ptycho2DTPA{ ii }.sol.rPIE_alpha;
%         metrics.rand_spos_subset_pct( ii ) = sim_ptycho2DTPA{ ii }.sol.spos.rand_spos_subset_pct;
        
    end
    
    name_data = num2str( sim_ptycho2DTPA{ii}.sol.rPIE_alpha, 'rPIE_alpha = %0.8f');
    
    metrics.meas_all_avg_of_log = meas_all_avg;
    metrics.name_data           = name_data;
    
    %=====================
    % create metrics plots
    %=====================

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

%     name_data = num2str( [ metrics.rPIE_alpha( 1 ), metrics.rand_spos_subset_pct( 1 ) ], 'rPIE_alpha = %0.8f, MBpct = %0.4f');
     
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
%     xlabel('Epoch')
%     ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
%     hold off
%     title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
%     grid on
%     ylim( y_lim )
%     legend( 'Location', 'southwest' ) 


end

%====================================================================================================================================================

































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
%     y_lim = [ -1, 6 ];
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
%     name_data = num2str( sim_ptycho2DTPA{ii}.sol.rPIE_alpha, 'rPIE_alpha = %0.8f');
%     
%     xlabel('Epoch')
%     ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
%     hold off
%     title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
%     grid on
%     ylim( y_lim )
%     legend( 'Location', 'southwest' ) 
% 
% 
% end
% 
% %====================================================================================================================================================



