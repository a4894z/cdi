%

clear; 

close all;

%========

cdi_full_ER     = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/fullbatch_total_vs_stoch/cdi_full_ER/sim_ptycho2DTPA.mat' );
cdi_full_mRAAR  = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/fullbatch_total_vs_stoch/cdi_full_mRAAR/sim_ptycho2DTPA.mat' );
cdi_full_RAAR   = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/fullbatch_total_vs_stoch/cdi_full_RAAR/sim_ptycho2DTPA.mat' );
cdi_stoch_ER    = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/fullbatch_total_vs_stoch/cdi_stoch_ER/sim_ptycho2DTPA.mat' );
cdi_stoch_mRAAR = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/fullbatch_total_vs_stoch/cdi_stoch_mRAAR/sim_ptycho2DTPA.mat' );
cdi_stoch_RAAR  = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/fullbatch_total_vs_stoch/cdi_stoch_RAAR/sim_ptycho2DTPA.mat' );

%========

skip = 5;

%========

figure; 

hold on
plot( cdi_full_mRAAR.sol.it.mtot( 1 : skip : end ), log10( cdi_full_mRAAR.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
plot( cdi_stoch_mRAAR.sol.it.mtot( 1 : skip : end ), log10( cdi_stoch_mRAAR.sol.metrics.meas_all( 1 : skip : end ) ), '--', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
plot( cdi_full_RAAR.sol.it.mtot( 1 : skip : end ), log10( cdi_full_RAAR.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.7 ] )
plot( cdi_stoch_RAAR.sol.it.mtot( 1 : skip : end ), log10( cdi_stoch_RAAR.sol.metrics.meas_all( 1 : skip : end ) ), '--', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.7 ] )
xlabel('Epoch')
ylabel('Cost Function Value')
hold off
title('$ log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s, ~ 100\% $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% legend({'cdi\_full\_ER', 'cdi\_full\_mRAAR', 'cdi\_full\_RAAR', 'cdi\_stoch\_ER', 'cdi\_stoch\_mRAAR', 'cdi\_stoch\_RAAR' })
legend({'cdi\_full\_mRAAR', 'cdi\_stoch\_mRAAR', 'cdi\_full\_RAAR', 'cdi\_stoch\_RAAR' })
grid on

%========

figure; 

hold on
plot( cdi_full_mRAAR.sol.it.mtot( 1 : skip : end ), log10( cdi_full_mRAAR.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
plot( cdi_stoch_mRAAR.sol.it.mtot( 1 : skip : end ), log10( cdi_stoch_mRAAR.sol.metrics.meas_all( 1 : skip : end ) ), '--', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
plot( cdi_full_ER.sol.it.mtot( 1 : skip : end ), log10( cdi_full_ER.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
plot( cdi_stoch_ER.sol.it.mtot( 1 : skip : end ), log10( cdi_stoch_ER.sol.metrics.meas_all( 1 : skip : end ) ), '--', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
xlabel('Epoch')
ylabel('Cost Function Value')
hold off
title('$ log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s, ~ 100\% $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% legend({'cdi\_full\_ER', 'cdi\_full\_mRAAR', 'cdi\_full\_RAAR', 'cdi\_stoch\_ER', 'cdi\_stoch\_mRAAR', 'cdi\_stoch\_RAAR' })
legend({'cdi\_full\_mRAAR', 'cdi\_stoch\_mRAAR', 'cdi\_full\_ER', 'cdi\_stoch\_ER' })
grid on

%========

% skip = 100;
% 
% figure; 
% 
% hold on
% plot( cdi_full_ER.sol.it.mtot( 1 : skip : end ), cdi_full_ER.sol.epoch_t( cdi_full_ER.sol.it.mtot( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
% plot( cdi_full_mRAAR.sol.it.mtot( 1 : skip : end ), cdi_full_mRAAR.sol.epoch_t( cdi_full_mRAAR.sol.it.mtot( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
% plot( cdi_full_RAAR.sol.it.mtot( 1 : skip : end ), cdi_full_RAAR.sol.epoch_t( cdi_full_RAAR.sol.it.mtot( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.7 ] )
% plot( cdi_stoch_ER.sol.it.mtot( 1 : skip : end ), cdi_stoch_ER.sol.epoch_t( cdi_stoch_ER.sol.it.mtot( 1 : skip : end ) ), '--', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.7 ] )
% plot( cdi_stoch_mRAAR.sol.it.mtot( 1 : skip : end ), cdi_stoch_mRAAR.sol.epoch_t( cdi_stoch_mRAAR.sol.it.mtot( 1 : skip : end ) ), '--', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.7 ] )
% plot( cdi_stoch_RAAR.sol.it.mtot( 1 : skip : end ), cdi_stoch_RAAR.sol.epoch_t( cdi_stoch_RAAR.sol.it.mtot( 1 : skip : end ) ), '--', 'linewidth', 2, 'color', [ 0.7, 0.5, 0.0 ] )
% xlabel('Epoch')
% ylabel('Time Elapsed Using Tic-Toc, Seconds')
% hold off
% legend({'cdi\_full\_ER', 'cdi\_full\_mRAAR', 'cdi\_full\_RAAR', 'cdi\_stoch\_ER', 'cdi\_stoch\_mRAAR', 'cdi\_stoch\_RAAR' })
% grid on

%====================================================================================================================================================

cdi_full_mRAAR_r2  = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/fullbatch_total_vs_stoch_r2/cdi_full_mRAAR_r2/sim_ptycho2DTPA.mat' );
cdi_stoch_mRAAR_r2 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/fullbatch_total_vs_stoch_r2/cdi_stoch_mRAAR_r2/sim_ptycho2DTPA.mat' );
cdi_full_RAAR_r2  = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/fullbatch_total_vs_stoch_r2/cdi_full_RAAR_r2/sim_ptycho2DTPA.mat' );
cdi_stoch_RAAR_r2 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/fullbatch_total_vs_stoch_r2/cdi_stoch_RAAR_r2/sim_ptycho2DTPA.mat' );
cdi_full_ER_r2  = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/fullbatch_total_vs_stoch_r2/cdi_full_ER_r2/sim_ptycho2DTPA.mat' );
cdi_stoch_ER_r2 = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/fullbatch_total_vs_stoch_r2/cdi_stoch_ER_r2/sim_ptycho2DTPA.mat' );

%========

skip = 5;

%========

figure; 

hold on
plot( cdi_full_mRAAR_r2.sol.it.mtot( 1 : skip : end ), log10( cdi_full_mRAAR_r2.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
plot( cdi_stoch_mRAAR_r2.sol.it.mtot( 1 : skip : end ), log10( cdi_stoch_mRAAR_r2.sol.metrics.meas_all( 1 : skip : end ) ), '--', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
plot( cdi_full_RAAR_r2.sol.it.mtot( 1 : skip : end ), log10( cdi_full_RAAR_r2.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.7 ] )
plot( cdi_stoch_RAAR_r2.sol.it.mtot( 1 : skip : end ), log10( cdi_stoch_RAAR_r2.sol.metrics.meas_all( 1 : skip : end ) ), '--', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.7 ] )
xlabel('Epoch')
ylabel('Cost Function Value')
hold off
title('$ log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s, ~ 100\% $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% legend({'cdi\_full\_ER', 'cdi\_full\_mRAAR', 'cdi\_full\_RAAR', 'cdi\_stoch\_ER', 'cdi\_stoch\_mRAAR', 'cdi\_stoch\_RAAR' })
legend({'cdi\_full\_mRAAR', 'cdi\_stoch\_mRAAR', 'cdi\_full\_RAAR', 'cdi\_stoch\_RAAR' })
grid on

%========

figure; 

hold on
plot( cdi_full_mRAAR_r2.sol.it.mtot( 1 : skip : end ), log10( cdi_full_mRAAR_r2.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
plot( cdi_stoch_mRAAR_r2.sol.it.mtot( 1 : skip : end ), log10( cdi_stoch_mRAAR_r2.sol.metrics.meas_all( 1 : skip : end ) ), '--', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
plot( cdi_full_ER_r2.sol.it.mtot( 1 : skip : end ), log10( cdi_full_ER_r2.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
plot( cdi_stoch_ER_r2.sol.it.mtot( 1 : skip : end ), log10( cdi_stoch_ER_r2.sol.metrics.meas_all( 1 : skip : end ) ), '--', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
xlabel('Epoch')
ylabel('Cost Function Value')
hold off
title('$ log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s, ~ 100\% $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% legend({'cdi\_full\_ER', 'cdi\_full\_mRAAR', 'cdi\_full\_RAAR', 'cdi\_stoch\_ER', 'cdi\_stoch\_mRAAR', 'cdi\_stoch\_RAAR' })
legend({'cdi\_full\_mRAAR', 'cdi\_stoch\_mRAAR', 'cdi\_full\_ER', 'cdi\_stoch\_ER' })
grid on

% %========
% 
% skip = 10;
% 
% figure; 
% 
% title('Exitwave Update')
% hold on
% plot( cdi_full_ER_r2.sol.it.mtot( 1 : skip : end ), cdi_full_ER_r2.sol.dt_exwv_up( cdi_full_ER_r2.sol.it.mtot( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
% plot( cdi_stoch_ER_r2.sol.it.mtot( 1 : skip : end ), cdi_stoch_ER_r2.sol.dt_exwv_up( cdi_stoch_ER_r2.sol.it.mtot( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
% plot( cdi_full_mRAAR_r2.sol.it.mtot( 1 : skip : end ), cdi_full_mRAAR_r2.sol.dt_exwv_up( cdi_full_mRAAR_r2.sol.it.mtot( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.7 ] )
% plot( cdi_stoch_mRAAR_r2.sol.it.mtot( 1 : skip : end ), cdi_stoch_mRAAR_r2.sol.dt_exwv_up( cdi_stoch_mRAAR_r2.sol.it.mtot( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.7 ] )
% xlabel('Epoch')
% ylabel('Time Elapsed Using Tic-Toc, Seconds')
% hold off
% legend({'cdi\_full\_ER', 'cdi\_stoch\_ER', 'cdi\_full\_mRAAR', 'cdi\_stoch\_mRAAR'})
% grid on
% 
% %========
% 
% skip = 10;
% 
% figure; 
% 
% title('Sample Update')
% hold on
% % plot( cdi_full_ER_r2.sol.it.mtot( 1 : skip : end ), cdi_full_ER_r2.sol.dt_sample_up( cdi_full_ER_r2.sol.it.mtot( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
% % plot( cdi_stoch_ER_r2.sol.it.mtot( 1 : skip : end ), cdi_stoch_ER_r2.sol.dt_sample_up( cdi_stoch_ER_r2.sol.it.mtot( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
% plot( cdi_full_mRAAR_r2.sol.it.mtot( 1 : skip : end ), cdi_full_mRAAR_r2.sol.dt_sample_up( cdi_full_mRAAR_r2.sol.it.mtot( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.7 ] )
% plot( cdi_stoch_mRAAR_r2.sol.it.mtot( 1 : skip : end ), cdi_stoch_mRAAR_r2.sol.dt_sample_up( cdi_stoch_mRAAR_r2.sol.it.mtot( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.7 ] )
% xlabel('Epoch')
% ylabel('Time Elapsed Using Tic-Toc, Seconds')
% hold off
% legend({ 'cdi\_full\_mRAAR', 'cdi\_stoch\_mRAAR'})
% grid on
% 
% %========
% 
% skip = 10;
% 
% figure; 
% 
% title('Probe Update')
% hold on
% % plot( cdi_full_ER_r2.sol.it.mtot( 1 : skip : end ), cdi_full_ER_r2.sol.dt_sample_up( cdi_full_ER_r2.sol.it.mtot( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.0 ] )
% % plot( cdi_stoch_ER_r2.sol.it.mtot( 1 : skip : end ), cdi_stoch_ER_r2.sol.dt_sample_up( cdi_stoch_ER_r2.sol.it.mtot( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
% plot( cdi_full_mRAAR_r2.sol.it.mtot( 1 : skip : end ), cdi_full_mRAAR_r2.sol.dt_probe_up( cdi_full_mRAAR_r2.sol.it.mtot( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.7 ] )
% plot( cdi_stoch_mRAAR_r2.sol.it.mtot( 1 : skip : end ), cdi_stoch_mRAAR_r2.sol.dt_probe_up( cdi_stoch_mRAAR_r2.sol.it.mtot( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.7 ] )
% xlabel('Epoch')
% ylabel('Time Elapsed Using Tic-Toc, Seconds')
% hold off
% legend({ 'cdi\_full\_mRAAR', 'cdi\_stoch\_mRAAR'})
% grid on



% plot( cdi_full_mRAAR_r2.sol.it.mtot( 1 : skip : end ), cdi_full_mRAAR_r2.sol.epoch_t( cdi_full_mRAAR_r2.sol.it.mtot( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
% plot( cdi_full_RAAR.sol.it.mtot( 1 : skip : end ), cdi_full_RAAR.sol.epoch_t( cdi_full_RAAR.sol.it.mtot( 1 : skip : end ) ), '-', 'linewidth', 2, 'color', [ 0.0, 0.0, 0.7 ] )
% plot( cdi_stoch_ER.sol.it.mtot( 1 : skip : end ), cdi_stoch_ER.sol.epoch_t( cdi_stoch_ER.sol.it.mtot( 1 : skip : end ) ), '--', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.7 ] )
% plot( cdi_stoch_mRAAR.sol.it.mtot( 1 : skip : end ), cdi_stoch_mRAAR.sol.epoch_t( cdi_stoch_mRAAR.sol.it.mtot( 1 : skip : end ) ), '--', 'linewidth', 2, 'color', [ 0.7, 0.0, 0.7 ] )
% plot( cdi_stoch_RAAR.sol.it.mtot( 1 : skip : end ), cdi_stoch_RAAR.sol.epoch_t( cdi_stoch_RAAR.sol.it.mtot( 1 : skip : end ) ), '--', 'linewidth', 2, 'color', [ 0.7, 0.5, 0.0 ] )


