%
%{

clear; close all; testing_mb05_rPIE_vs_alpha

%}

%====================================================================================================================================================

addpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi/misc/output/shadedErrorBar/' );  

%====================================================================================================================================================

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/mb0p05/';

%====================================================================================================================================================



























%====================================================================================================================================================

path_data = {};
N_trials  = [];

%========

% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p0001/independenttrials_08Sep2021_t235455/' ];
% N_trials( end + 1 )  = 8;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p001/independenttrials_08Sep2021_t121600/' ];
% N_trials( end + 1 )  = 10;
% 
% % path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p0025/independenttrials_10Sep2021_t092114/' ];
% % N_trials( end + 1 )  = 10;
% % 
% % path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p005/independenttrials_09Sep2021_t214215/' ];
% % N_trials( end + 1 )  = 10;
% % 
% % path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p0075/independenttrials_09Sep2021_t095311/' ];
% % N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p01/independenttrials_16Aug2021_t173355/' ];
% N_trials( end + 1 )  = 10;
% 
% % path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p03/independenttrials_03Sep2021_t163043/' ];
% % N_trials( end + 1 )  = 10;
% % 
% % path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p05/independenttrials_17Aug2021_t105554/' ];
% % N_trials( end + 1 )  = 10;
% % 
% % path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p07/independenttrials_04Sep2021_t094428/' ];
% % N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p10/independenttrials_18Aug2021_t053649/' ];
% N_trials( end + 1 )  = 10;
% 
% % path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p30/independenttrials_19Aug2021_t013134/' ];
% % N_trials( end + 1 )  = 10;
% % 
% % path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p50/independenttrials_05Sep2021_t030345/' ];
% % N_trials( end + 1 )  = 10;
% % 
% % path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p70/independenttrials_05Sep2021_t154235/' ];
% % N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha1p00/independenttrials_06Sep2021_t033324/' ];
% N_trials( end + 1 )  = 10;

%========

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p001_oldsample_for_probe_update/independenttrials_04Nov2021_t042443/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p01_recompute_exitwaves_after_sample/independenttrials_24Sep2021_t230058/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p01_oldsample_for_probe_update/independenttrials_24Sep2021_t123331/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p10_oldsample_for_probe_update/independenttrials_23Oct2021_t210230/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha1p00_oldsample_for_probe_update/independenttrials_30Oct2021_t051441/' ];
N_trials( end + 1 )  = 10;


%====================================================================================================================================================

for jj = 1 : length( path_data )
    
    load_and_plot( path_data{ jj }, N_trials( jj ) )
    
end

%====================================================================================================================================================

function load_and_plot( path_data, N_trials )

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

    %========
    
    skip = 1;

    h1 = figure();  
    set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )

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
    
end

%====================================================================================================================================================







