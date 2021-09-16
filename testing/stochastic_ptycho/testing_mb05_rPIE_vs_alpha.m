%
%{

clear; close all; testing_mb05_rPIE_vs_alpha

%}

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/mb0p05/';

path_data = {};
N_trials  = [];

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p0001/independenttrials_08Sep2021_t235455/' ];
N_trials( end + 1 )  = 8;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p001/independenttrials_08Sep2021_t121600/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p0025/independenttrials_10Sep2021_t092114/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p005/independenttrials_09Sep2021_t214215/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p0075/independenttrials_09Sep2021_t095311/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p01/independenttrials_16Aug2021_t173355/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p03/independenttrials_03Sep2021_t163043/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p05/independenttrials_17Aug2021_t105554/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p07/independenttrials_04Sep2021_t094428/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p10/independenttrials_18Aug2021_t053649/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p30/independenttrials_19Aug2021_t013134/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p50/independenttrials_05Sep2021_t030345/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p70/independenttrials_05Sep2021_t154235/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha1p00/independenttrials_06Sep2021_t033324/' ];
N_trials( end + 1 )  = 10;

%========

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

        meas_all_avg = meas_all_avg + log10( sim_ptycho2DTPA{ ii }.sol.metrics.meas_all ); 

    end

    meas_all_avg = meas_all_avg / N_trials;

    %========

    skip = 1;

%     figure; 
    h1 = figure();  
    set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )

    hold on

    for ii = 1 : N_trials

        plot( sim_ptycho2DTPA{ ii }.sol.it.mtot( 1 : skip : end ), log10( sim_ptycho2DTPA{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

    end

    plot( sim_ptycho2DTPA{ ii }.sol.it.mtot( 1 : skip : end ), meas_all_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )

    name_data = num2str( [ sim_ptycho2DTPA{ii}.sol.rPIE_alpha, sim_ptycho2DTPA{ii}.sol.spos.rand_spos_subset_pct ], 'rPIE_alpha = %0.8f, MBpct = %0.4f');
    
    xlabel('Epoch')
    ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
    hold off
    title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
    grid on
    ylim( y_lim )
    legend( 'Location', 'northeast' ) 


end

%====================================================================================================================================================







