
close all; clear;

%{

restoredefaultpath
addpath( genpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi' ));   



clear; close all; testing_sample_projection_phaseramping

%}



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