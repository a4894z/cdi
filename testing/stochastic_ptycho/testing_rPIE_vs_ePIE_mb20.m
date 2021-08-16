%
%{

clear; close all; testing_rPIE_vs_ePIE_mb20

%}

%====================================================================================================================================================

% alpha = 0.01

%==========================================================
% These are WITH sample magnitude less than one constraints
%==========================================================

for ii = 1 : 10
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/03Aug2021_withTleq1/cdi_ER_rPIE_0p01_0p99_gpu4/independenttrials_03Aug2021_t105121/' ;
    cdi_ER_rPIE_0p01_0p99{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : 10
    
    plot( cdi_ER_rPIE_0p01_0p99{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p01_0p99{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('rPIE ( 0.01, 0.99 ), Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s,~ 20\% $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim([0, 6])
legend 


%=============================================================
% These are WITHOUT sample magnitude less than one constraints
%=============================================================

for ii = 1 : 10
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/02Aug2021_withoutTleq1/cdi_ER_rPIE_0p01_0p99_gpu1/independenttrials_02Aug2021_t113401/' ;
    cdi_ER_rPIE_0p01_0p99{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : 10
    
    plot( cdi_ER_rPIE_0p01_0p99{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p01_0p99{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('rPIE ( 0.01, 0.99 ), Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s,~ 20\% $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim([0, 6])


%====================================================================================================================================================

% alpha = 0.05

%========

for ii = 1 : 10
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/02Aug2021_withoutTleq1/cdi_ER_rPIE_0p05_0p95_gpu3/independenttrials_02Aug2021_t110532/' ;
    cdi_ER_rPIE_0p05_0p95{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : 10
    
    plot( cdi_ER_rPIE_0p05_0p95{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p05_0p95{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('rPIE ( 0.05, 0.95 ), Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s,~ 20\% $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim([0, 6])

%=============================================================
% These are WITHOUT sample magnitude less than one constraints
%=============================================================

for ii = 1 : 10
    
    path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/03Aug2021_withTleq1/cdi_ER_rPIE_0p05_0p95_gpu3/independenttrials_03Aug2021_t105052/' ;
    cdi_ER_rPIE_0p05_0p95{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

end

skip = 1;

figure; 

hold on

for ii = 1 : 10
    
    plot( cdi_ER_rPIE_0p05_0p95{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p05_0p95{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

end

xlabel('Epoch')
ylabel('rPIE ( 0.05, 0.95 ), Cost Function Value')
hold off
title('$log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s,~ 20\% $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
grid on
ylim([0, 6])

%========






























% %=======================================================
% % These are WITH probe photons constraint and shrinkwrap
% %=======================================================
% 
% for ii = 1 : 10
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/02Aug2021_withoutTleq1/cdi_ER_ePIE_gpu4/independenttrials_01Aug2021_t165028/' ;
%     cdi_ER_ePIE{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : 10
%     
%     plot( cdi_ER_ePIE{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_ePIE{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% xlabel('Epoch')
% ylabel('ePIE, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s,~ 20\% $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim([0, 6])
% 
% %========
% 
% for ii = 1 : 10
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/02Aug2021_withoutTleq1/cdi_ER_rPIE_0p75_0p25_gpu2/independenttrials_01Aug2021_t164914/' ;
%     cdi_ER_rPIE_0p75_0p25{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : 10
%     
%     plot( cdi_ER_rPIE_0p75_0p25{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p75_0p25{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% xlabel('Epoch')
% ylabel('rPIE ( 0.75, 0.25 ), Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s,~ 20\% $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim([0, 6])
% 
% %========
% 
% for ii = 1 : 10
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/02Aug2021_withoutTleq1/cdi_ER_rPIE_0p5_0p5_gpu3/independenttrials_01Aug2021_t164958/' ;
%     cdi_ER_rPIE_0p5_0p5{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : 10
%     
%     plot( cdi_ER_rPIE_0p5_0p5{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p5_0p5{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% xlabel('Epoch')
% ylabel('rPIE ( 0.5, 0.5 ), Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s,~ 20\% $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim([0, 6])
% 
% %========
% 
% for ii = 1 : 10
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/02Aug2021_withoutTleq1/cdi_ER_rPIE_0p25_0p75_gpu2/independenttrials_02Aug2021_t105805/' ;
%     cdi_ER_rPIE_0p25_0p75{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : 10
%     
%     plot( cdi_ER_rPIE_0p25_0p75{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p25_0p75{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% xlabel('Epoch')
% ylabel('rPIE ( 0.25, 0.75 ), Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s,~ 20\% $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim([0, 6])
% 
% %========
% 
% for ii = 1 : 10
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/02Aug2021_withoutTleq1/cdi_ER_rPIE_0p15_0p85_gpu4/independenttrials_02Aug2021_t110220/' ;
%     cdi_ER_rPIE_0p15_0p85{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : 10
%     
%     plot( cdi_ER_rPIE_0p15_0p85{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p15_0p85{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% xlabel('Epoch')
% ylabel('rPIE ( 0.15, 0.85 ), Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg] ~\forall ~N_s,~ 20\% $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim([0, 6])
% 







%==========================================================
% These are WITHOUT probe photons constraint and shrinkwrap
%==========================================================












%====================================================================================================================================================


%========


%========
