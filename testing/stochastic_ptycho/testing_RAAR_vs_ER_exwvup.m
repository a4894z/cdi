% 

clear; close all;

%====================================================================================================================================================

cdi_RAAR_01 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/RAAR_vs_ER_Exwvup/cdi_RAAR_epoch500_exwv1/sim_ptycho2DTPA.mat' );
cdi_RAAR_05 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/RAAR_vs_ER_Exwvup/cdi_RAAR_epoch500_exwv5/sim_ptycho2DTPA.mat' );
cdi_RAAR_10 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/RAAR_vs_ER_Exwvup/cdi_RAAR_epoch500_exwv10/sim_ptycho2DTPA.mat' );
% cdi_RAAR_20 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/RAAR_vs_ER_Exwvup/cdi_RAAR_epoch500_exwv20/sim_ptycho2DTPA.mat' );
% cdi_RAAR_50 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/RAAR_vs_ER_Exwvup/cdi_RAAR_epoch500_exwv50/sim_ptycho2DTPA.mat' );

%========

figure; 

hold on
semilogy( cdi_RAAR_01.sol.it.mtot, log10( cdi_RAAR_01.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0, 0, 0 ] )
semilogy( cdi_RAAR_05.sol.it.mtot, log10( cdi_RAAR_05.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.8, 0, 0 ] )
semilogy( cdi_RAAR_10.sol.it.mtot, log10( cdi_RAAR_10.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0, 0.8, 0 ] )
% semilogy( cdi_RAAR_20.sol.it.mtot, log10( cdi_RAAR_20.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0, 0, 0.8 ] )
% semilogy( cdi_RAAR_50.sol.it.mtot, log10( cdi_RAAR_50.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0, 0.8, 0.8 ] )
hold off
title('$ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F ~\forall ~N_s $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
legend({'1', ...
        '5', ...
        '10' })

%====================================================================================================================================================

cdi_ER_01 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/RAAR_vs_ER_Exwvup/cdi_ER_epoch500_exwv1/sim_ptycho2DTPA.mat' );
cdi_ER_05 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/RAAR_vs_ER_Exwvup/cdi_ER_epoch500_exwv5/sim_ptycho2DTPA.mat' );

figure; 

hold on
semilogy( cdi_ER_01.sol.it.mtot, log10( cdi_ER_01.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0, 0, 0 ] )
semilogy( cdi_ER_05.sol.it.mtot, log10( cdi_ER_05.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.8, 0, 0 ] )
% semilogy( cdi_ER_10.sol.it.mtot, log10( cdi_ER_10.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0, 0.8, 0 ] )
% semilogy( cdi_ER_20.sol.it.mtot, log10( cdi_ER_20.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0, 0, 0.8 ] )
% semilogy( cdi_ER_50.sol.it.mtot, log10( cdi_ER_50.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0, 0.8, 0.8 ] )
hold off
title('$ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F ~\forall ~N_s $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
legend({'1', ...
        '5', ...
        '10' })
