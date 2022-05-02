

clear; close all; clc;

% load /net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/processed_expt/L0258_to_L0265_combined_256x768.mat;

new_mat = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/processed_expt/L0274_to_L0280_combined_512x512.mat' );
old_mat = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/processed_expt/old/L0274_to_L0280_combined_512x512.mat' );


                       
                       

norm( old_mat.expt.spos.xh_yv_zb_frameL( : ) - new_mat.expt.spos.xh_yv_zb_frameL(:) )

norm( old_mat.expt.spos.yv_xh_frameL( : ) - new_mat.expt.spos.yv_xh_frameL(:) )






tmp0 = new_mat.expt.spos.rsALL;

s_x = -1.15;
sol.spos.shear_x = [ [ 1,   0 ]; ...
                     [ s_x, 1 ] ];

tmp0 = transpose( sol.spos.shear_x * transpose( tmp0 ));

tmp0( :, 1 ) = tmp0( :, 1 ) + 1;
norm( old_mat.expt.spos.rsALL( : ) - tmp0(:) )



figure; plot(old_mat.expt.spos.rsALL(:,1))
figure; plot(1 + tmp0(:,1))



figure; plot(tmp0(:,1) - old_mat.expt.spos.rsALL(:,1))


close all

figure; 
plot_2Dscan_positions( tmp0, [], [], [] )
set( gca, 'xdir', 'reverse' )
set( gca, 'ydir', 'normal' )
xlabel('xh, lab frame'); ylabel('yv, lab frame');
daspect([1 1 1])  

figure; 
plot_2Dscan_positions( old_mat.expt.spos.rsALL, [], [], [] )
set( gca, 'xdir', 'reverse' )
set( gca, 'ydir', 'normal' )
xlabel('xh, lab frame'); ylabel('yv, lab frame');
daspect([1 1 1])  














tmp0 = old_mat.expt.spos.rs;

s_x = 1.15;
sol.spos.shear_x = [ [ 1,   0 ]; ...
                     [ s_x, 1 ] ];

tmp0 = transpose( sol.spos.shear_x * transpose( tmp0 ));



close all

figure; plot( old_mat.expt.spos.rs )
figure; plot( tmp0 )






figure; plot( expt.spos.rs )



















tmp0 = expt.spos.yv_xh_frameL( expt.spos.indxsubset, : );
tmp1 = old_mat.expt.spos.yv_xh_frameL( old_mat.expt.spos.indxsubset, : );


expt.spos.yv_xh_frameL( :, 2 ) = -1 * expt.spos.yv_xh_frameL( :, 2 );
tmp0( :, 2 ) =  -1 * tmp0( :, 2 );

close all

figure; 
plot_2Dscan_positions( expt.spos.yv_xh_frameL, [], tmp0( 1 : 1 : end, : ), [] )
set( gca, 'xdir', 'reverse' )
set( gca, 'ydir', 'normal' )
xlabel('xh, lab frame'); ylabel('yv, lab frame');
daspect([1 1 1])  

figure; 
plot_2Dscan_positions( old_mat.expt.spos.yv_xh_frameL, [], tmp1( 1 : 1 : end, : ), [] )
set( gca, 'xdir', 'reverse' )
set( gca, 'ydir', 'normal' )
xlabel('xh, lab frame'); ylabel('yv, lab frame');
daspect([1 1 1])  
