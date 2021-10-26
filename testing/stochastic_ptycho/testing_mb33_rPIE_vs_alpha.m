% 

%
%{

clear; close all; testing_mb33_rPIE_vs_alpha

%}

%====================================================================================================================================================

addpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/shadedErrorBar/' );  

%====================================================================================================================================================

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/mb0p33/';

%====================================================================================================================================================

path_data = {};
N_trials  = [];

%========

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha1p00_oldsample_for_probe_update/independenttrials_07Oct2021_t204113/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p80_oldsample_for_probe_update/independenttrials_07Oct2021_t002735/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p60_oldsample_for_probe_update/independenttrials_06Oct2021_t040251/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p40_oldsample_for_probe_update/independenttrials_05Oct2021_t075837/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p20_oldsample_for_probe_update/independenttrials_04Oct2021_t115500/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p10_oldsample_for_probe_update/independenttrials_03Oct2021_t155218/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p02_oldsample_for_probe_update/independenttrials_22Oct2021_t163410/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p01_oldsample_for_probe_update/independenttrials_25Sep2021_t132343/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p001_oldsample_for_probe_update/independenttrials_26Sep2021_t002144/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p0001_oldsample_for_probe_update/independenttrials_26Sep2021_t112108/' ];
N_trials( end + 1 )  = 10;

%====================================================================================================================================================

for jj = 1 : length( path_data )
    
    metrics{ jj } = load_and_plot( path_data{ jj }, N_trials( jj ) );   %#ok<SAGROW>
    
    metrics{ jj }.rootpath_data = rootpath_data;                        %#ok<SAGROW>

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

%     name_data = num2str( [ sim_ptycho2DTPA{ 1 }.sol.rPIE_alpha, sim_ptycho2DTPA{ 1 }.sol.spos.rand_spos_subset_pct ], 'rPIE_alpha = %0.8f, MBpct = %0.4f');

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

%====================================================================================================================================================





























% %========
% 
% for jj = 1 : length( path_data )
%     
%     load_and_plot( path_data{ jj }, N_trials( jj ) )
%     
%     5;
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
% end
% 
% %====================================================================================================================================================







