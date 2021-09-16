% 

%
%{

clear; close all; testing_mb50_rPIE_vs_alpha

%}


rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/mb0p50/';

path_data = {};
N_trials  = [];

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

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p5_alpha0p05/independenttrials_06Sep2021_t062131/' ];
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











% 
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb20_alpha100/independenttrials_11Aug2021_t083622/';
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
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb20_alpha80/independenttrials_15Aug2021_t234638/';
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
% clear('cdi_ER_rPIE_0p80')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb20_alpha60/independenttrials_15Aug2021_t061743/';
%     cdi_ER_rPIE_0p60{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p60_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p60_avg = cdi_ER_rPIE_0p60_avg + log10( cdi_ER_rPIE_0p60{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p60_avg = cdi_ER_rPIE_0p60_avg / N_trials;
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
%     plot( cdi_ER_rPIE_0p60{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p60{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p60{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p60_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.60, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear('cdi_ER_rPIE_0p60')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb20_alpha40/independenttrials_14Aug2021_t130526/';
%     cdi_ER_rPIE_0p40{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p40_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p40_avg = cdi_ER_rPIE_0p40_avg + log10( cdi_ER_rPIE_0p40{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p40_avg = cdi_ER_rPIE_0p40_avg / N_trials;
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
%     plot( cdi_ER_rPIE_0p40{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p40{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p40{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p40_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.40, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear('cdi_ER_rPIE_0p40')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb20_alpha20/independenttrials_13Aug2021_t202125/';
%     cdi_ER_rPIE_0p20{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p20_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p20_avg = cdi_ER_rPIE_0p20_avg + log10( cdi_ER_rPIE_0p20{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p20_avg = cdi_ER_rPIE_0p20_avg / N_trials;
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
%     plot( cdi_ER_rPIE_0p20{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p20{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p20{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p20_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.20, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'northeast' ) 
% 
% %========
% 
% clear('cdi_ER_rPIE_0p20')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb20_alpha10/independenttrials_13Aug2021_t084624/';
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
% clear('cdi_ER_rPIE_0p10')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb20_alpha09/independenttrials_12Aug2021_t224717/';
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
% clear('cdi_ER_rPIE_0p09')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb20_alpha08/independenttrials_12Aug2021_t130610/';
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
% clear('cdi_ER_rPIE_0p08')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb20_alpha07/independenttrials_12Aug2021_t033840/';
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
% clear('cdi_ER_rPIE_0p07')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb20_alpha06/independenttrials_11Aug2021_t181332/';
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
% clear('cdi_ER_rPIE_0p06')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb20_alpha05_gpu3/independenttrials_10Aug2021_t071308/';
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
% clear('cdi_ER_rPIE_0p05')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb20_alpha04_gpu3/independenttrials_09Aug2021_t214747/';
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
% clear('cdi_ER_rPIE_0p04')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb20_alpha03_gpu3/independenttrials_09Aug2021_t081023/';
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
% clear('cdi_ER_rPIE_0p03')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb20_alpha02_gpu3/independenttrials_06Aug2021_t175923/';
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
% clear('cdi_ER_rPIE_0p02')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb20_alpha01_gpu2/independenttrials_04Aug2021_t145912/';
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
% clear('cdi_ER_rPIE_0p01')
% 
% %====================================================================================================================================================