%
%{

clear; close all; testing_mb_alpha0p01

%}

y_lim = [-1, 6];

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb01_alpha01/independenttrials_24Aug2021_t083135/';
    cdi_ER_rPIE_mb01_0p01{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

%========

cdi_ER_rPIE_mb01_0p01_avg = 0;

for ii = 1 : N_trials
    
    cdi_ER_rPIE_mb01_0p01_avg = cdi_ER_rPIE_mb01_0p01_avg + cdi_ER_rPIE_mb01_0p01{ ii }.sol.metrics.meas_all;

end

cdi_ER_rPIE_mb01_0p01_avg = cdi_ER_rPIE_mb01_0p01_avg / N_trials;

%========

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_rPIE_mb01_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_mb01_0p01{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

plot( cdi_ER_rPIE_mb01_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_mb01_0p01_avg( 1 : skip : end ) ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )

xlabel('Epoch')
ylabel('mb0.01, rPIE alpha = 0.01, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'northeast' ) 

%========

clear( 'cdi_ER_rPIE_mb01_0p01' )

%====================================================================================================================================================

N_trials = 9;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/gpu3/cdi_rPIE_mb005_alpha01/independenttrials_24Aug2021_t083220/';
    cdi_ER_rPIE_mb005_0p01{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

%========

cdi_ER_rPIE_mb005_0p01_avg = 0;

for ii = 1 : N_trials
    
    cdi_ER_rPIE_mb005_0p01_avg = cdi_ER_rPIE_mb005_0p01_avg + cdi_ER_rPIE_mb005_0p01{ ii }.sol.metrics.meas_all;

end

cdi_ER_rPIE_mb005_0p01_avg = cdi_ER_rPIE_mb005_0p01_avg / N_trials;

%========

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_rPIE_mb005_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_mb005_0p01{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

plot( cdi_ER_rPIE_mb005_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_mb005_0p01_avg( 1 : skip : end ) ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )

xlabel('Epoch')
ylabel('mb0.005, rPIE alpha = 0.01, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'southwest' ) 

%========

clear( 'cdi_ER_rPIE_mb005_0p01' )

%====================================================================================================================================================

N_trials = 10;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb05_alpha01_gpu2/independenttrials_16Aug2021_t173355/';
    cdi_ER_rPIE_mb05_0p01{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

%========

cdi_ER_rPIE_mb05_0p01_avg = 0;

for ii = 1 : N_trials
    
    cdi_ER_rPIE_mb05_0p01_avg = cdi_ER_rPIE_mb05_0p01_avg + cdi_ER_rPIE_mb05_0p01{ ii }.sol.metrics.meas_all;

end

cdi_ER_rPIE_mb05_0p01_avg = cdi_ER_rPIE_mb05_0p01_avg / N_trials;

%========

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_rPIE_mb05_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_mb05_0p01{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

plot( cdi_ER_rPIE_mb05_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_mb05_0p01_avg( 1 : skip : end ) ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )

xlabel('Epoch')
ylabel('mb0.05, rPIE alpha = 0.01, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'northeast' ) 

%========

clear( 'cdi_ER_rPIE_mb05_0p01' )

%====================================================================================================================================================

N_trials = 8;

for ii = 1 : N_trials
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_mb025_alpha01/independenttrials_24Aug2021_t235115/';
    cdi_ER_rPIE_mb025_0p01{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

%========

cdi_ER_rPIE_mb025_0p01_avg = 0;

for ii = 1 : N_trials
    
    cdi_ER_rPIE_mb025_0p01_avg = cdi_ER_rPIE_mb025_0p01_avg + cdi_ER_rPIE_mb025_0p01{ ii }.sol.metrics.meas_all;

end

cdi_ER_rPIE_mb025_0p01_avg = cdi_ER_rPIE_mb025_0p01_avg / N_trials;

%========

skip = 1;

figure; 

hold on

for ii = 1 : N_trials
    
    plot( cdi_ER_rPIE_mb025_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_mb025_0p01{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

plot( cdi_ER_rPIE_mb025_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_mb025_0p01_avg( 1 : skip : end ) ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )

xlabel('Epoch')
ylabel('mb0.025, rPIE alpha = 0.01, Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim( y_lim )
legend( 'Location', 'northeast' ) 

%========

clear( 'cdi_ER_rPIE_mb025_0p01' )

%====================================================================================================================================================

