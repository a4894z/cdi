% 

%
%{

clear; close all; testing_mb10_rPIE_vs_alpha

%}

y_lim = [-1, 6];

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha01_gpu3/independenttrials_06Aug2021_t095015/';
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
legend( 'Location', 'northeast' ) 

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha02_gpu2/independenttrials_09Aug2021_t164930/';
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
legend( 'Location', 'northeast' ) 

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha03_gpu2/independenttrials_10Aug2021_t031506/';
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
legend( 'Location', 'northeast' ) 

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha04_gpu2/independenttrials_10Aug2021_t113245/';
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
legend( 'Location', 'northeast' ) 

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha05_gpu2/independenttrials_10Aug2021_t195103/';
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
legend( 'Location', 'northeast' ) 

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha06_gpu2/independenttrials_11Aug2021_t035951/';
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
legend( 'Location', 'northeast' ) 

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha07_gpu2/independenttrials_11Aug2021_t121230/';
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
legend( 'Location', 'northeast' ) 

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha08_gpu2/independenttrials_11Aug2021_t202809/';
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
legend( 'Location', 'northeast' ) 

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha09_gpu2/independenttrials_12Aug2021_t044127/';
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
legend( 'Location', 'northeast' ) 

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha10_gpu2/independenttrials_12Aug2021_t125535/';
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
legend( 'Location', 'northeast' ) 

%====================================================================================================================================================














