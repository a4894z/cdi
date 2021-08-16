% 

%
%{

close all; clear; testing_stochGD_rPIE_vs_alpha

%}

y_lim = [-1, 6];

%====================================================================================================================================================

for ii = 1 : 10
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_stochGD_alpha01_gpu3/independenttrials_04Aug2021_t150026/' ;
    cdi_ER_rPIE_0p01{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : 10
    
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

for ii = 1 : 10
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_stochGD_alpha02_gpu4/independenttrials_05Aug2021_t095936/' ;
    cdi_ER_rPIE_0p02{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : 10
    
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

for ii = 1 : 10
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_stochGD_alpha03_gpu2/independenttrials_05Aug2021_t100113/' ;
    cdi_ER_rPIE_0p03{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : 10
    
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


for ii = 1 : 10
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_stochGD_alpha04_gpu3/independenttrials_05Aug2021_t105850/' ;
    cdi_ER_rPIE_0p04{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : 10
    
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

for ii = 1 : 10
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_stochGD_alpha05_gpu4/independenttrials_06Aug2021_t093611/' ;
    cdi_ER_rPIE_0p05{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : 10
    
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

for ii = 1 : 10
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_stochGD_alpha06/independenttrials_09Aug2021_t214746/' ;
    cdi_ER_rPIE_0p06{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : 10
    
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

for ii = 1 : 10
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_stochGD_alpha07/independenttrials_10Aug2021_t130027/' ;
    cdi_ER_rPIE_0p07{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : 10
    
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

for ii = 1 : 10
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_stochGD_alpha08/independenttrials_11Aug2021_t042043/' ;
    cdi_ER_rPIE_0p08{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : 10
    
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

for ii = 1 : 10
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_stochGD_alpha09/independenttrials_11Aug2021_t194747/' ;
    cdi_ER_rPIE_0p09{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : 10
    
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


for ii = 1 : 10
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_stochGD_alpha10/independenttrials_12Aug2021_t111313/' ;
    cdi_ER_rPIE_0p10{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : 10
    
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


