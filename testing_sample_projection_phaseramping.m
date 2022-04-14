

%{

cd '/net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi';

restoredefaultpath
addpath( genpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi' ));   



clear; close all; testing_sample_projection_phaseramping

%}




data_mat_path = '/net/s8iddata/export/8-id-ECA/pmyint/2022_Analysis/Multislice_Simulations/Laminography_simulations/PetraIII_Eiffel_simulation_lami_1deg/';



kk = 1;

for ss = 0 : 20 : ( 360 - 1 )
    
    tmp0 = load( [ data_mat_path, num2str( ss, 'Trans_proj%.5d.mat' ) ]);
    
    
    sz = size( tmp0.data );
    slice_sz_r = 64;
    slice_sz_c = 1024;

    localslice_r = (( 0.5 * sz( 1 ) + 1 ) - 0.5 * slice_sz_r ) : (( 0.5 * sz( 1 ) + 0 ) + 0.5 * slice_sz_r );
    localslice_c = (( 0.5 * sz( 2 ) + 1 ) - 0.5 * slice_sz_c ) : (( 0.5 * sz( 2 ) + 0 ) + 0.5 * slice_sz_c );

    localslice_r = round( localslice_r );
    localslice_c = round( localslice_c );

    
    figure( 666 ); 
    
    ax1 = subplot(121);
    imagesc( abs( tmp0.data( localslice_r, localslice_c ) ), [ 0.8, 1.0 ] )
%     daspect([1 1 1])
    colormap( ax1, 'winter' )
    colorbar
    
    ax2 = subplot(122);
    imagesc( angle( tmp0.data( localslice_r, localslice_c ) ), 0.7 * [ -pi, pi ] )
%     daspect([1 1 1])
    colormap( ax2, 'hsv' )
    colorbar 
%     subplot(133)
%     imagescHSV( tmp0.data )
%     daspect([1 1 1])
    
    title( num2str( ss, 'rotation angle index = %d' ))
    drawnow
    
%     subplot(122)
%     plot( sampleTF_2Dproj_localview( 1 : end, 772, ss ) )
    
    tmp1 = tmp0.data( localslice_r, localslice_c );
    the_sinogram( :, kk ) = tmp1( : );
    kk = kk + 1;

end



return





data_hdf_path = '/net/s8iddata/export/8-id-ECA/pmyint/2022_Analysis/Multislice_Simulations/Laminography_simulations/PetraIII_Eiffel_simulation_lami_1deg/';
data_hdf_name = 'PETRAIII_Eiffel_transf_funproj_lami_1degree.h5';

data_hdf = [ data_hdf_path, data_hdf_name ];

h5info = h5info( data_hdf );

sampleTF_2Dproj = h5read( [ data_hdf_path, data_hdf_name ], '/exchange/data' );
sampleTF_2Dproj = permute( sampleTF_2Dproj, [ 2, 1, 3 ] );

sz = size( sampleTF_2Dproj );
slice_sz_r = 768 + 64;
slice_sz_c = 768 + 64;

localslice_r = (( 0.5 * sz( 1 ) + 1 ) - 0.5 * slice_sz_r ) : (( 0.5 * sz( 1 ) + 0 ) + 0.5 * slice_sz_r );
localslice_c = (( 0.5 * sz( 2 ) + 1 ) - 0.5 * slice_sz_c ) : (( 0.5 * sz( 2 ) + 0 ) + 0.5 * slice_sz_c );

localslice_r = round( localslice_r );
localslice_c = round( localslice_c );

sampleTF_2Dproj_localview = sampleTF_2Dproj( localslice_r, localslice_c, : );

for ss = 80 %1 : 10 : sz( 3 )
    
    figure( 666 ); 
%     subplot(121)
    imagesc( sampleTF_2Dproj_localview( :, :, ss ), [ 0, 1 ] )
    daspect([1 1 1])
    colormap winter
    title( num2str( [ ss, sampleTF_2Dproj_localview( 607, 185, ss ) ], 'rotation angle index = %d, substrate value = %.4f' ))
    drawnow
    
%     subplot(122)
%     plot( sampleTF_2Dproj_localview( 1 : end, 772, ss ) )
    
    

end


return


tmp0 = load( '/net/s8iddata/export/8-id-ECA/pmyint/2022_Analysis/2022_phase_beating _study/transfer_fun_projections/Thickness95nm_angle1p01_goldonsilicon.mat' );
sample_proj{ 1 } = tmp0.projection;

tmp0 = load( '/net/s8iddata/export/8-id-ECA/pmyint/2022_Analysis/2022_phase_beating _study/transfer_fun_projections/Thickness95nm_angle1p11_goldonsilicon.mat' );
sample_proj{ 2 }  = tmp0.projection;

tmp0 = load( '/net/s8iddata/export/8-id-ECA/pmyint/2022_Analysis/2022_phase_beating _study/transfer_fun_projections/Thickness95nm_angle1p21_goldonsilicon.mat' );
sample_proj{ 3 } = tmp0.projection;

tmp0 = load( '/net/s8iddata/export/8-id-ECA/pmyint/2022_Analysis/2022_phase_beating _study/transfer_fun_projections/Thickness95nm_angle1p31_goldonsilicon.mat' );
sample_proj{ 4 } = tmp0.projection;

tmp0 = load( '/net/s8iddata/export/8-id-ECA/pmyint/2022_Analysis/2022_phase_beating _study/transfer_fun_projections/Thickness95nm_angle1p41_goldonsilicon.mat' );
sample_proj{ 5 } = tmp0.projection;

tmp0 = load( '/net/s8iddata/export/8-id-ECA/pmyint/2022_Analysis/2022_phase_beating _study/transfer_fun_projections/Thickness95nm_angle1p51_goldonsilicon.mat' );
sample_proj{ 6 } = tmp0.projection;

tmp0 = load( '/net/s8iddata/export/8-id-ECA/pmyint/2022_Analysis/2022_phase_beating _study/transfer_fun_projections/Thickness95nm_angle1p61_goldonsilicon.mat' );
sample_proj{ 7 } = tmp0.projection;

tmp0 = load( '/net/s8iddata/export/8-id-ECA/pmyint/2022_Analysis/2022_phase_beating _study/transfer_fun_projections/Thickness95nm_angle1p71_goldonsilicon.mat' );
sample_proj{ 8 } = tmp0.projection;

tmp0 = load( '/net/s8iddata/export/8-id-ECA/pmyint/2022_Analysis/2022_phase_beating _study/transfer_fun_projections/Thickness95nm_angle1p81_goldonsilicon.mat' );
sample_proj{ 9 } = tmp0.projection;

tmp0 = load( '/net/s8iddata/export/8-id-ECA/pmyint/2022_Analysis/2022_phase_beating _study/transfer_fun_projections/Thickness95nm_angle1p91_goldonsilicon.mat' );
sample_proj{ 10 } = tmp0.projection;


sz0 = size( sample_proj{ 1 } );



for ii = 1 : length( sample_proj )
    
    sz = size( sample_proj{ ii } );
    
    slcr = round( ( 0.5 * sz( 1 ) - 0.5 * sz0( 1 ) + 1 ) : ( 0.5 * sz( 1 ) + 0.5 * sz0( 1 ) ));
    slcc = round( ( 0.5 * sz( 2 ) - 0.5 * sz0( 2 ) + 1 ) : ( 0.5 * sz( 2 ) + 0.5 * sz0( 2 ) ));



    h1 = figure(); 
    set( h1, 'Visible', 'on', 'Position',[ 10, 1, 700, 1080 ] )
    
    tmp0 = sample_proj{ ii };
    imagescHSV( tmp0( slcr, slcc ) )
    daspect([1 1 1])
    
    export_fig( num2str( ii, 'proj_%d' ), '-r120.0' )
    close all;
    
    
end