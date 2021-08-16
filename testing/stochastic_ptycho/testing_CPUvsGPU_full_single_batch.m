%

clear; close all;

%====================================================================================================================================================

cdi_full_CPU = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/r2_fullbatch_stochcoord/cdi_fullbatch_CPU/sim_ptycho2DTPA.mat' );
cdi_full_GPU = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/r2_fullbatch_stochcoord/cdi_fullbatch_GPU/sim_ptycho2DTPA.mat' );

cdi_single_CPU = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/r2_fullbatch_stochcoord/cdi_stochcoord_CPU/sim_ptycho2DTPA.mat' );
cdi_single_GPU = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/r2_fullbatch_stochcoord/cdi_stochcoord_GPU/sim_ptycho2DTPA.mat' );

%========

close all;

figure; 

hold on
semilogy( cdi_full_CPU.sol.it.mtot, log10( cdi_full_CPU.sol.metrics.meas_all ), '-', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.0 ] )
semilogy( cdi_full_GPU.sol.it.mtot, log10( cdi_full_GPU.sol.metrics.meas_all ), '-', 'linewidth', 2, 'color', [ 0.8, 0.0, 0.0 ] )
semilogy( cdi_single_CPU.sol.it.mtot, log10( cdi_single_CPU.sol.metrics.meas_all ), '-', 'linewidth', 2, 'color', [ 0.0, 0.8, 0.0 ] )
semilogy( cdi_single_GPU.sol.it.mtot, log10( cdi_single_GPU.sol.metrics.meas_all ), '-', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.8 ] )
hold off
grid on
title('$ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F ~\forall ~N_s $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
legend({'CPU Full','GPU Full','CPU Stoch','GPU Stoch'})
xlabel('Epoch')
ylabel('Cost')

figure; 

subplot(211)
hold on
semilogy( cdi_full_CPU.sol.epoch_t, '-', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.0 ] )
semilogy( cdi_full_GPU.sol.epoch_t, '-', 'linewidth', 2, 'color', [ 0.8, 0.0, 0.0 ] )
semilogy( cdi_single_CPU.sol.epoch_t, '-', 'linewidth', 2, 'color', [ 0.0, 0.8, 0.0 ] )
semilogy( cdi_single_GPU.sol.epoch_t, '-', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.8 ] )
hold off
grid on
legend({'CPU Full','GPU Full','CPU Stoch','GPU Stoch'})
xlabel('Epoch')
ylabel('Time [s]')

subplot(212)
hold on
semilogy( cdi_full_CPU.sol.epoch_t ./ cdi_full_GPU.sol.epoch_t, '-', 'linewidth', 2, 'color', [ 0.8, 0.0, 0.0 ] )
semilogy( cdi_single_CPU.sol.epoch_t ./ cdi_single_GPU.sol.epoch_t, '-', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.8 ] )
hold off
legend({'CPU / GPU Full','CPU / GPU Stoch'})
xlabel('Epoch')
ylabel('CPU/GPU Speedup')
grid on







%========
% 
% figure; 
% 
% hold on
% semilogy( cdi_full_CPU.sol.it.mtot, log10( cdi_full_CPU.sol.metrics.meas_all ), '-', 'linewidth', 2, 'color', [ 0, 0, 0 ] )
% semilogy( cdi_full_GPU.sol.it.mtot, log10( cdi_full_GPU.sol.metrics.meas_all ), '-', 'linewidth', 2, 'color', [ 0.8, 0, 0 ] )
% hold off
% grid on
% title('$ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F ~\forall ~N_s $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% legend({'CPU','GPU'})
% xlabel('Epoch')
% ylabel('Cost')
% 
% figure; 
% 
% subplot(211)
% hold on
% semilogy( cdi_full_CPU.sol.epoch_t, '-', 'linewidth', 2, 'color', [ 0, 0, 0 ] )
% semilogy( cdi_full_GPU.sol.epoch_t, '-', 'linewidth', 2, 'color', [ 0.8, 0, 0.8 ] )
% hold off
% grid on
% legend({'CPU','GPU'})
% xlabel('Epoch')
% ylabel('Time [s]')
% 
% subplot(212)
% semilogy( cdi_full_CPU.sol.epoch_t ./ cdi_full_GPU.sol.epoch_t, '-', 'linewidth', 2, 'color', [ 0.0, 0.8, 0.0 ] )
% xlabel('Epoch')
% ylabel('CPU/GPU Speedup')
% grid on
% 
% %========
% 
% figure; 
% 
% hold on
% semilogy( cdi_single_CPU.sol.it.mtot, log10( cdi_single_CPU.sol.metrics.meas_all ), '-', 'linewidth', 2, 'color', [ 0, 0, 0 ] )
% semilogy( cdi_single_GPU.sol.it.mtot, log10( cdi_single_GPU.sol.metrics.meas_all ), '-', 'linewidth', 2, 'color', [ 0.8, 0, 0 ] )
% hold off
% grid on
% title('$ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F ~\forall ~N_s $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% legend({'CPU','GPU'})
% xlabel('Epoch')
% ylabel('Cost')
% 
% figure; 
% 
% subplot(211)
% hold on
% semilogy( cdi_single_CPU.sol.epoch_t, '-', 'linewidth', 2, 'color', [ 0, 0, 0 ] )
% semilogy( cdi_single_GPU.sol.epoch_t, '-', 'linewidth', 2, 'color', [ 0.8, 0, 0.8 ] )
% hold off
% grid on
% legend({'CPU','GPU'})
% xlabel('Epoch')
% ylabel('Time [s]')
% 
% subplot(212)
% semilogy( cdi_single_CPU.sol.epoch_t ./ cdi_single_GPU.sol.epoch_t, '-', 'linewidth', 2, 'color', [ 0.0, 0.8, 0.0 ] )
% xlabel('Epoch')
% ylabel('CPU/GPU Speedup')
% grid on
% 












%====================================================================================================================================================










% 
% 
% 
% cdi_full_CPU   = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/r1_fullbatch_stochcoord/cdi_fullbatch_CPU/sim_ptycho2DTPA.mat' );
% cdi_full_GPU   = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/r1_fullbatch_stochcoord/cdi_fullbatch_GPU/sim_ptycho2DTPA.mat' );
% cdi_single_CPU = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/r1_fullbatch_stochcoord/cdi_stochcoord_CPU/sim_ptycho2DTPA.mat' );
% cdi_single_GPU = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/r1_fullbatch_stochcoord/cdi_stochcoord_GPU/sim_ptycho2DTPA.mat' );
% 
% %========
% 
% figure; 
% 
% hold on
% semilogy( cdi_full_CPU.sol.it.mtot, log10( cdi_full_CPU.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0, 0, 0 ] )
% semilogy( cdi_full_GPU.sol.it.mtot, log10( cdi_full_GPU.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.8, 0, 0 ] )
% hold off
% title('$ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F ~\forall ~N_s $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% 
% 
% figure; 
% 
% subplot(211)
% hold on
% semilogy( cdi_full_CPU.sol.epoch_t, '-o', 'linewidth', 2, 'color', [ 0, 0, 0 ] )
% semilogy( cdi_full_GPU.sol.epoch_t, '-o', 'linewidth', 2, 'color', [ 0.8, 0, 0.8 ] )
% hold off
% 
% subplot(212)
% semilogy( cdi_full_CPU.sol.epoch_t ./ cdi_full_GPU.sol.epoch_t, '-o', 'linewidth', 2, 'color', [ 0.8, 0, 0.8 ] )
% 
% %========
% 
% figure; 
% 
% hold on
% semilogy( cdi_single_CPU.sol.it.mtot, log10( cdi_single_CPU.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0, 0, 0 ] )
% semilogy( cdi_single_GPU.sol.it.mtot, log10( cdi_single_GPU.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.8, 0, 0 ] )
% hold off
% title('$ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F ~\forall ~N_s $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% xlabel('Epoch')
% ylabel('Cost')
% 
% 
% figure; 
% 
% subplot(211)
% hold on
% semilogy( cdi_single_CPU.sol.epoch_t, '-o', 'linewidth', 2, 'color', [ 0, 0, 0 ] )
% semilogy( cdi_single_GPU.sol.epoch_t, '-o', 'linewidth', 2, 'color', [ 0.8, 0, 0.8 ] )
% hold off
% legend({'CPU','GPU'})
% xlabel('Epoch')
% ylabel('Time [s]')
% 
% subplot(212)
% semilogy( cdi_single_CPU.sol.epoch_t ./ cdi_single_GPU.sol.epoch_t, '-o', 'linewidth', 2, 'color', [ 0.8, 0, 0.8 ] )
% xlabel('Epoch')
% ylabel('CPU/GPU Speedup')
