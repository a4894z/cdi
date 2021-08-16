
clear; close all;

cdi_40 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/fixed_epoch_vs_minibatch/r2/cdi_40/sim_ptycho2DTPA.mat' );
cdi_20 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/fixed_epoch_vs_minibatch/r2/cdi_20/sim_ptycho2DTPA.mat' );
cdi_10 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/fixed_epoch_vs_minibatch/r2/cdi_10/sim_ptycho2DTPA.mat' );

cdi_5  = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/fixed_epoch_vs_minibatch/r2/cdi_5/sim_ptycho2DTPA.mat' );
cdi_4  = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/fixed_epoch_vs_minibatch/r2/cdi_4/sim_ptycho2DTPA.mat' );
cdi_3  = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/fixed_epoch_vs_minibatch/r2/cdi_3/sim_ptycho2DTPA.mat' );
cdi_2  = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/fixed_epoch_vs_minibatch/r2/cdi_2/sim_ptycho2DTPA.mat' );
cdi_1a = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/fixed_epoch_vs_minibatch/r2/cdi_1a/sim_ptycho2DTPA.mat' );
cdi_1  = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/fixed_epoch_vs_minibatch/r2/cdi_1/sim_ptycho2DTPA.mat' );

%========

figure; 

hold on
plot( cdi_5.sol.it.mtot, log10( cdi_5.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
plot( cdi_4.sol.it.mtot, log10( cdi_4.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
plot( cdi_3.sol.it.mtot, log10( cdi_3.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.7 ] )
plot( cdi_2.sol.it.mtot, log10( cdi_2.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.7, 0.7, 0.0 ] )
plot( cdi_1.sol.it.mtot, log10( cdi_1.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.7 ] )
hold off
title('$ log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
legend({'5%', '4%', '3%', '2%', '1%'})
grid on

%========

figure; 

hold on
plot( cdi_40.sol.it.mtot, log10( cdi_40.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
plot( cdi_20.sol.it.mtot, log10( cdi_20.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
plot( cdi_10.sol.it.mtot, log10( cdi_10.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.7 ] )
plot( cdi_5.sol.it.mtot, log10( cdi_5.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.7, 0.7, 0.0 ] )
plot( cdi_1.sol.it.mtot, log10( cdi_1.sol.metrics.meas_all ), '-o', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.7 ] )
hold off
title('$ log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
legend({'40%', '20%', '10%', '5%', '1%'})
grid on

%========

tmp0 = cdi_20;

figure; 
hold on

for bb = 1 : size( tmp0.sol.metrics.meas, 2 )
    
    plot( tmp0.sol.it.mtot, log10( tmp0.sol.metrics.meas( :, bb ) ), '-o', 'linewidth', 2 )

end

hold off
title('$ log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
legend
grid on
xlabel('Epoch')
ylabel('Cost')

%========

tmp0 = cdi_5;

figure; 
hold on

for bb = 1 : size( tmp0.sol.metrics.meas, 2 )
    
    plot( tmp0.sol.it.mtot, log10( tmp0.sol.metrics.meas( :, bb ) ), '-o', 'linewidth', 2 )

end

hold off
title('$ log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
legend
grid on
xlabel('Epoch')
ylabel('Cost')

%====================================================================================================================================================

figure; 

hold on
plot( cdi_5.sol.epoch_t , '-o', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
plot( cdi_4.sol.epoch_t , '-o', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
plot( cdi_3.sol.epoch_t , '-o', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.7 ] )
plot( cdi_2.sol.epoch_t , '-o', 'linewidth', 2, 'color', [ 0.7, 0.7, 0.0 ] )
plot( cdi_1.sol.epoch_t , '-o', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.7 ] )
hold off
xlabel('Epoch')
ylabel('Time')
legend({'5%', '4%', '3%', '2%', '1%'})
grid on

%========

figure; 

hold on
plot( cdi_40.sol.epoch_t , '-o', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
plot( cdi_20.sol.epoch_t , '-o', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
plot( cdi_10.sol.epoch_t , '-o', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.7 ] )
plot( cdi_5.sol.epoch_t , '-o', 'linewidth', 2, 'color', [ 0.7, 0.7, 0.0 ] )
plot( cdi_1.sol.epoch_t , '-o', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.7 ] )
hold off
xlabel('Epoch')
ylabel('Time')
legend({'40%', '20%', '10%', '5%', '1%'})
grid on

%====================================================================================================================================================

cdi_10_CPU = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/CPU_vs_GPU/cdi_10_CPU/sim_ptycho2DTPA.mat' );
cdi_10_GPU = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/CPU_vs_GPU/cdi_10_GPU/sim_ptycho2DTPA.mat' );

figure; 

hold on
plot( cdi_10_CPU.sol.epoch_t , '-o', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
plot( cdi_10_GPU.sol.epoch_t , '-o', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
hold off
xlabel('Epoch')
ylabel('Time')
% legend({'40%', '20%', '10%', '5%', '1%'})
grid on

figure; 

hold on
plot( cdi_10_CPU.sol.epoch_t ./ cdi_10_GPU.sol.epoch_t, '-o', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
hold off
xlabel('Epoch')
ylabel('CPU / GPU Timing Ratio')
legend({'10%'})
% legend({'40%', '20%', '10%', '5%', '1%'})
grid on

%====================================================================================================================================================

cdi_100_CPU = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/CPU_vs_GPU/cdi_fullbatch_CPU/sim_ptycho2DTPA.mat' );
cdi_100_GPU = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/stochastic_ptycho_benchmarking/CPU_vs_GPU/cdi_fullbatch_GPU/sim_ptycho2DTPA.mat' );

figure; 

hold on
plot( cdi_100_CPU.sol.epoch_t , '-o', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
plot( cdi_100_GPU.sol.epoch_t , '-o', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
hold off
xlabel('Epoch')
ylabel('Time')
% legend({'40%', '20%', '10%', '5%', '1%'})
grid on

figure; 

hold on
plot( cdi_100_CPU.sol.epoch_t ./ cdi_100_GPU.sol.epoch_t, '-o', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
hold off
xlabel('Epoch')
ylabel('CPU / GPU Timing Ratio')
legend({'100%'})
% legend({'40%', '20%', '10%', '5%', '1%'})
grid on




figure; 

hold on
plot( cdi_100_CPU.sol.epoch_t ./ cdi_100_GPU.sol.epoch_t, '-o', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.0 ] )
plot( cdi_10_CPU.sol.epoch_t ./ cdi_10_GPU.sol.epoch_t, '-o', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )

hold off
xlabel('Epoch')
ylabel('CPU / GPU Timing Ratio')
legend({'100%', '10%'})
% legend({'40%', '20%', '10%', '5%', '1%'})
grid on










