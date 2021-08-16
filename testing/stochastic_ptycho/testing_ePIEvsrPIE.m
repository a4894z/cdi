
%==================================================================================================

clear; close all;
load /home/ash/Documents/Science/Matlab/Code/cdi/1-run_no_probe_fro2_or_occ_constraint/G0024_apslogo_27May2020_t132403_it25001.mat;

nm = 3; 
aalpha = 0.8; 

tmp0 = abs( probe.P( :, :, nm )).^2; 
tmp1 = max( tmp0(:) );

u_ePIE = -tmp0 + (1/aalpha) * tmp1; 

figure; imagesc( u_ePIE ); colormap jet; colorbar; title('ePIE'); daspect([1 1 1])


u_rPIE = -aalpha * tmp0 + aalpha * tmp1; 

figure; imagesc( u_rPIE ); colormap jet; colorbar; title('rPIE'); daspect([1 1 1])

norm( u_ePIE, 'fro' )
norm( u_rPIE, 'fro' )


%==================================================================================================