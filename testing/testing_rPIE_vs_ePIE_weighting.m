
%
%{

restoredefaultpath; 
addpath( genpath( pwd ));
clearvars -except expt sol




clear; close all; testing_rPIE_vs_ePIE_weighting


%}

%====================================================================================================================================================

load /net/s8iddata/export/8-id-ECA/Analysis/atripath/gpu4/cdi/L0258_to_L0265_combined_512x512.mat

%====================================================================================================================================================

probe = sol.probe.phi;

probe_intensity     = sum( abs( probe ) .^ 2, 3 );
max_probe_intensity = max( probe_intensity( : ) );

%====================================================================================================================================================

figure;

%========

aalpha = 1e-0;
w_rPIE = probe_intensity ./ ( ( 1 - aalpha ) * probe_intensity + aalpha * max_probe_intensity );

subplot(141);
imagesc( w_rPIE );
daspect([1 1 1]);
colorbar
colormap turbo
title( num2str( aalpha, 'r/ePIE weighting, alpha = %.1e' ))

%========

aalpha = 1e-3;
w_rPIE = probe_intensity ./ ( ( 1 - aalpha ) * probe_intensity + aalpha * max_probe_intensity );

subplot(142);
imagesc( w_rPIE );
daspect([1 1 1]);
colorbar
colormap turbo
title( num2str( aalpha, 'rPIE weighting, alpha = %.1e' ))

%========

aalpha = 1e-6;
w_rPIE = probe_intensity ./ ( ( 1 - aalpha ) * probe_intensity + aalpha * max_probe_intensity );

subplot(143);
imagesc( w_rPIE );
daspect([1 1 1]);
colorbar
colormap turbo
title( num2str( aalpha, 'rPIE weighting, alpha = %.1e' ))

%========

aalpha = 1e-9;
w_rPIE = probe_intensity ./ ( ( 1 - aalpha ) * probe_intensity + aalpha * max_probe_intensity );

subplot(144);
imagesc( w_rPIE );
daspect([1 1 1]);
colorbar
colormap turbo
title( num2str( aalpha, 'rPIE weighting, alpha = %.1e' ))

%====================================================================================================================================================

figure;

%========

aalpha = 1e-0;

u_ePIE = ( 1 / aalpha ) * max_probe_intensity - probe_intensity;
u_rPIE = aalpha * ( max_probe_intensity - probe_intensity );

subplot(231);
imagesc( u_rPIE );
daspect([1 1 1]);
colorbar
colormap turbo
title( num2str( aalpha, 'rPIE penalty, alpha = %.1e' ))

subplot(234);
imagesc( u_ePIE );
daspect([1 1 1]);
colorbar
colormap turbo
title( num2str( aalpha, 'ePIE penalty, alpha = %.1e' ))

%========

aalpha = 1e-3;

u_ePIE = ( 1 / aalpha ) * max_probe_intensity - probe_intensity;
u_rPIE = aalpha * ( max_probe_intensity - probe_intensity );

subplot(232);
imagesc( u_rPIE );
daspect([1 1 1]);
colorbar
colormap turbo
title( num2str( aalpha, 'rPIE penalty, alpha = %.1e' ))

subplot(235);
imagesc( u_ePIE );
daspect([1 1 1]);
colorbar
colormap turbo
title( num2str( aalpha, 'ePIE penalty, alpha = %.1e' ))

%========

aalpha = 1e-6;

u_ePIE = ( 1 / aalpha ) * max_probe_intensity - probe_intensity;
u_rPIE = aalpha * ( max_probe_intensity - probe_intensity );

subplot(233);
imagesc( u_rPIE );
daspect([1 1 1]);
colorbar
colormap turbo
title( num2str( aalpha, 'rPIE penalty, alpha = %.1e' ))

subplot(236);
imagesc( u_ePIE );
daspect([1 1 1]);
colorbar
colormap turbo
title( num2str( aalpha, 'ePIE penalty, alpha = %.1e' ))

%========

