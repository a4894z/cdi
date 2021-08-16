%

cdi_orthmodes_1 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/modes_orthog_vary/cdi_probe_orthog_study_1/sim_ptycho2DTPA.mat' );
cdi_orthmodes_10  = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/modes_orthog_vary/cdi_probe_orthog_study_10/sim_ptycho2DTPA.mat' );
% cdi_orthmodes_50 = load( '' );
cdi_orthmodes_100 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/modes_orthog_vary/cdi_probe_orthog_study_100/sim_ptycho2DTPA.mat' );

%========

figure; 

hold on
plot( cdi_orthmodes_1.sol.it.mtot, log10( cdi_orthmodes_1.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
plot( cdi_orthmodes_10.sol.it.mtot, log10( cdi_orthmodes_10.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
plot( cdi_orthmodes_100.sol.it.mtot, log10( cdi_orthmodes_100.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.7 ] )
xlabel('Epoch')
ylabel('Cost Function Value')
hold off
title('$ log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s, ~ 10\% $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
legend({'1', '10', '100'})
grid on
