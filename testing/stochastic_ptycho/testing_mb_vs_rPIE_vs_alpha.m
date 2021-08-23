
%
%{

clear; close all; testing_mb_vs_rPIE_vs_alpha

%}

y_lim = [-1, 6];

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_stochGD_alpha01_gpu3/independenttrials_04Aug2021_t150026/';
    cdi_stoch_rPIE_0p01{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_stoch_rPIE_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_stoch_rPIE_0p01{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('stoch, rPIE alpha = 0.01, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'southwest' ) 

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb05_alpha01_gpu2/independenttrials_16Aug2021_t173355/';
    cdi_mb05_rPIE_0p01{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_mb05_rPIE_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_mb05_rPIE_0p01{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('mb05, rPIE alpha = 0.01, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'northeast' ) 

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb10_alpha01_gpu2/independenttrials_13Aug2021_t154650/';
    cdi_mb10_rPIE_0p01{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_mb10_rPIE_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_mb10_rPIE_0p01{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('mb10, rPIE alpha = 0.01, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'northeast' ) 

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb20_alpha01_gpu2/independenttrials_04Aug2021_t145912/';
    cdi_mb20_rPIE_0p01{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_mb20_rPIE_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_mb20_rPIE_0p01{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('mb20, rPIE alpha = 0.01, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'northeast' ) 

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb33_alpha01_gpu4/independenttrials_14Aug2021_t090209/';
    cdi_mb33_rPIE_0p01{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_mb33_rPIE_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_mb33_rPIE_0p01{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('mb33, rPIE alpha = 0.01, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'northeast' ) 

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha01_randT/independenttrials_14Aug2021_t081939/';
    cdi_mb100_rPIE_0p01{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_mb100_rPIE_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_mb100_rPIE_0p01{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('full, rPIE alpha = 0.01, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'northeast' ) 




