% 

%
%{

clear; close all; testing_fullGD_rPIE_vs_alpha

%}

y_lim = [-1, 6];

%======================================================================================================================
% using random sample starts, full batch update order doesn't matter so we need to introduce stochasticity in some form
%======================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_ePIE_full_randT/independenttrials_13Aug2021_t155848/';
    cdi_ER_ePIE{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_ePIE{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_ePIE{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('ePIE, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'southwest' )

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha01_randT/independenttrials_14Aug2021_t081939/';
    cdi_ER_rPIE_0p01{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_rPIE_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p01{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('rPIE alpha = 0.01, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'southwest' )

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha05_randT/independenttrials_15Aug2021_t004151/';
    cdi_ER_rPIE_0p05{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_rPIE_0p05{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p05{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('rPIE alpha = 0.05, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'southwest' )

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha10_randT/independenttrials_15Aug2021_t170501/';
    cdi_ER_rPIE_0p10{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_rPIE_0p10{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p10{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('rPIE alpha = 0.10, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'southwest' )












return

%=======================================================
% just one trial for same random starts for sample/probe
%=======================================================

N_trials = 1;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha01_gpu4/independenttrials_04Aug2021_t150149/';
    cdi_ER_rPIE_0p01{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_rPIE_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p01{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('rPIE alpha = 0.01, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'southwest' )

%====================================================================================================================================================

N_trials = 1;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha02_gpu4/independenttrials_07Aug2021_t011459/';
    cdi_ER_rPIE_0p02{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_rPIE_0p02{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p02{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('rPIE alpha = 0.02, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'southwest' )

%====================================================================================================================================================

N_trials = 1;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha03_gpu4/independenttrials_07Aug2021_t174121/';
    cdi_ER_rPIE_0p03{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_rPIE_0p03{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p03{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('rPIE alpha = 0.03, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'southwest' )

%====================================================================================================================================================

N_trials = 1;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha04_gpu4/independenttrials_08Aug2021_t100521/';
    cdi_ER_rPIE_0p04{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_rPIE_0p04{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p04{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('rPIE alpha = 0.04, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'southwest' )

%====================================================================================================================================================

N_trials = 1;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha05_gpu4/independenttrials_09Aug2021_t023601/';
    cdi_ER_rPIE_0p05{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_rPIE_0p05{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p05{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('rPIE alpha = 0.05, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'southwest' )

%====================================================================================================================================================

N_trials = 1;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha06/independenttrials_09Aug2021_t133912/';
    cdi_ER_rPIE_0p06{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_rPIE_0p06{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p06{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('rPIE alpha = 0.06, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'southwest' )

%====================================================================================================================================================

N_trials = 1;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha07/independenttrials_09Aug2021_t151739/';
    cdi_ER_rPIE_0p07{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_rPIE_0p07{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p07{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('rPIE alpha = 0.07, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'southwest' )

%====================================================================================================================================================

N_trials = 1;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha08/independenttrials_09Aug2021_t165553/';
    cdi_ER_rPIE_0p08{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_rPIE_0p08{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p08{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('rPIE alpha = 0.08, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'southwest' )

%====================================================================================================================================================

N_trials = 1;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha09/independenttrials_09Aug2021_t183320/';
    cdi_ER_rPIE_0p09{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_rPIE_0p09{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p09{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('rPIE alpha = 0.09, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'southwest' )

%====================================================================================================================================================

N_trials = 1;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha10/independenttrials_09Aug2021_t201045/';
    cdi_ER_rPIE_0p10{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_rPIE_0p10{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p10{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('rPIE alpha = 0.10, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'southwest' )

%====================================================================================================================================================



