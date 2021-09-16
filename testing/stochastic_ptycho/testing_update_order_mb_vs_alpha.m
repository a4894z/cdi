%% 

%
%{

clear; close all; testing_update_order_mb_vs_alpha

%}

path_data = {};
N_trials  = [];
ylabel_data = {};

%========

% rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/mb0p20/';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p05/independenttrials_10Aug2021_t071308/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p05_oldsample_for_probe_update/independenttrials_15Sep2021_t120331/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: only Tk';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p05_recompute_exitwaves_after_sample/independenttrials_15Sep2021_t120733/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/mb0p10/';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p05/independenttrials_10Aug2021_t195103/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p05_oldsample_for_probe_update/independenttrials_15Sep2021_t114753/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p05_recompute_exitwaves_after_sample/independenttrials_15Sep2021_t193044/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========

for jj = 1 : length( path_data )
    
    load_and_plot( path_data{ jj }, N_trials( jj ), ylabel_data{ jj } )
    
    5;
    
end

%====================================================================================================================================================

function load_and_plot( path_data, N_trials, ylabel_data )

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
    
    name_data = [ ylabel_data, ', ', name_data ];
    
    
    xlabel('Epoch')
    ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
    hold off
    title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
    grid on
    ylim( y_lim )
    legend( 'Location', 'northeast' ) 

end

%====================================================================================================================================================


