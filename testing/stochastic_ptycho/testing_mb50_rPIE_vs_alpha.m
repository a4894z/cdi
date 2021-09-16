% 

%
%{

clear; close all; testing_mb50_rPIE_vs_alpha

%}

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/mb0p50/';

path_data = {};
N_trials  = [];

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p0000001/independenttrials_07Sep2021_t090731/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p000001/independenttrials_06Sep2021_t160435/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p00001/independenttrials_05Sep2021_t225627/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p0001/independenttrials_05Sep2021_t054812/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p001/independenttrials_03Sep2021_t171323/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p01/independenttrials_04Sep2021_t163059/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p03/independenttrials_05Sep2021_t131034/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p04/independenttrials_09Sep2021_t045047/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p05/independenttrials_06Sep2021_t062131/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p06/independenttrials_09Sep2021_t215940/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p07/independenttrials_06Sep2021_t233002/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p08/independenttrials_10Sep2021_t150917/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p10/independenttrials_07Sep2021_t163956/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p20/independenttrials_11Sep2021_t081609/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p40/independenttrials_12Sep2021_t011547/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p60/independenttrials_12Sep2021_t181604/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p80/independenttrials_13Sep2021_t111918/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha1p00/independenttrials_14Sep2021_t042415/' ];
N_trials( end + 1 )  = 10;

%========

for jj = 1 : length( path_data )
    
    load_and_plot( path_data{ jj }, N_trials( jj ) )
    
    5;
    
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


