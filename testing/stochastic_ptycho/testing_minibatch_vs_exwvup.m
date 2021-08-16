%
clear; close all;

cdi_e1000_exwv_1 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/minibatch10/cdi_e1000_exwv1/sim_ptycho2DTPA.mat' );
cdi_e200_exwv_5  = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/minibatch10/cdi_e200_exwv5/sim_ptycho2DTPA.mat' );
cdi_e100_exwv_10 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/minibatch10/cdi_e100_exwv10/sim_ptycho2DTPA.mat' );



%========

figure; 

hold on
plot( cdi_e1000_exwv_1.sol.it.mtot, log10( cdi_e1000_exwv_1.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
plot( cdi_e200_exwv_5.sol.it.mtot, log10( cdi_e200_exwv_5.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
plot( cdi_e100_exwv_10.sol.it.mtot, log10( cdi_e100_exwv_10.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.7 ] )
xlabel('Epoch')
ylabel('Cost Function Value')
hold off
title('$ log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s, ~ 10\% $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
legend({'Epoch=1000, psi up=1', 'Epoch=200, psi up=5', 'Epoch=100, psi up=10'})
grid on

%========

figure; 

hold on
plot( cdi_e1000_exwv_1.sol.it.mtot, cdi_e1000_exwv_1.sol.epoch_t( cdi_e1000_exwv_1.sol.it.mtot ), '-o', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
plot( cdi_e200_exwv_5.sol.it.mtot, cdi_e200_exwv_5.sol.epoch_t( cdi_e200_exwv_5.sol.it.mtot ), '-o', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
plot( cdi_e100_exwv_10.sol.it.mtot, cdi_e100_exwv_10.sol.epoch_t( cdi_e100_exwv_10.sol.it.mtot ), '-o', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.7 ] )
xlabel('Epoch')
ylabel('Time Elapsed Using Tic-Toc, Seconds')
hold off
legend({'Epoch=1000, psi up=1, Ttotal = 0.7847 hr', 'Epoch=200, psi up=5, Ttotal = 0.3148 hr', 'Epoch=100, psi up=10, , Ttotal = 0.2519 hr'})
grid on