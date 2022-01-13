%
%{

%========

restoredefaultpath; 

% addpath( genpath( '~/Documents/Science/Matlab/Code/cdi/' ));
addpath( genpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi' ));

addpath( genpath( pwd ));

%========

clear; close all; testing_multislice_actual_vs_virtual_probe

%========

%}

%====================================================================================================================================================
% load sample volume, etc for multislice forward problem
%=======================================================

load /net/s8iddata/export/8-id-ECA/Analysis/atripath/Petra_aprime_0p4sampit.mat

%================================
% coordinate system length scales
%================================

resolution = resolution * 1e-6;   % um to m

expt.csys.z2.dLy = single( resolution( 2 ));
expt.csys.z2.dLx = single( resolution( 3 ));
expt.csys.z2.dLz = single( resolution( 1 ));

%========

expt.energy = energy;                           % x-ray energy ( keV )    
expt.lambda = ( 12.4 / expt.energy ) * 1e-10;   % wavelength ( in meters )

%====================================================================================================================================================

% figure; 
% subplot(131); 
% imagesc( abs( squeeze( sum( data, 1 ))));
% title('sum along rows projection')
% xlabel( 'pages, x' )
% ylabel( 'columns, y' )
% daspect( [ 1, 1, 1 ] )
% 
% subplot(132); 
% imagesc( abs( squeeze( sum( data, 2 ))));
% title('sum along columns projection')
% xlabel( 'pages, x' )
% ylabel( 'rows, z' )
% daspect( [ 1, 1, 1 ] )
% 
% subplot(133); 
% imagesc( abs( squeeze( sum( data, 3 ))));
% title('sum along pages projection')
% xlabel( 'columns, y' )
% ylabel( 'rows, z' )
% daspect( [ 1, 1, 1 ] )

%============================================================================================
% change coordinate system so that rows are vertical, columns are horizontal, pages are depth
%============================================================================================

%             data1 = data - min( data( : ));
%             
%             data1 = permute( data1, [ 3, 2, 1 ]);
%             data1 = permute( data1, [ 2, 1, 3 ]);
% %             data1 = flip( data1, 1 );
%             
%             V = 100 * imresize3( data1, 0.25 );
%             
%             figure; 
%             params.isoslvl = 1e-6;
%             params.alphalvl = 0.35;
%             isosurface_phase_surface_nopermute( abs( V ), params ); % colormap gray
%             daspect([1 1 1])
%             axis tight
% 
% 
%             figure; 
%             params.isoslvl = 1e-6;
%             params.alphalvl = 0.35;
%             isosurface_phase_surface( abs( V ), params ); % colormap gray
%             daspect([1 1 1])
%             axis tight
%             
%             
%             
%             
% 
% data0 = padarray( data, [ 0, 1, 0 ], 'post' );
% 
% data1 = permute( data0, [ 3, 2, 1 ]);
% data2 = permute( data1, [ 2, 1, 3 ]);
% 
% expt.sample.indxrefr_data_raw = permute( expt.sample.indxrefr_data_raw, [ 3, 2, 1 ]);
% 
% % expt.sample.indxrefr_data_raw = permute( data, [ 1, 3, 2 ]);
% % expt.sample.indxrefr_data_raw = permute( expt.sample.indxrefr_data_raw, [ 3, 2, 1 ]);

%========

expt.sample.indxrefr_data_raw = permute( data, [ 3, 2, 1 ]);
expt.sample.indxrefr_data_raw = permute( expt.sample.indxrefr_data_raw, [ 2, 1, 3 ]);
expt.sample.indxrefr_data_raw = flip( expt.sample.indxrefr_data_raw, 1 );
            
% % odd length arrays suck, zero pad to be even
% expt.sample.indxrefr_data_raw = padarray( expt.sample.indxrefr_data_raw, [ 1, 1, 0 ], 0, 'post' );

%=================
% get problem size
%=================

expt.sz.sz = single( size( expt.sample.indxrefr_data_raw ));
expt.sz.r  = expt.sz.sz( 1 );
expt.sz.c  = expt.sz.sz( 2 );
expt.sz.p  = expt.sz.sz( 3 );

%========

figure; 

subplot(231); 
imagesc( imag( squeeze( sum( expt.sample.indxrefr_data_raw, 1 ) / expt.sz.r )));
set( gca, 'ydir', 'normal' )
title('sum along rows projection, imag part')
xlabel( 'pages, z' )
ylabel( 'columns, x' )
daspect( [ 1, 1, 1 ] )
colorbar

subplot(232); 
imagesc( imag( squeeze( sum( expt.sample.indxrefr_data_raw, 2 ) / expt.sz.c )));
set( gca, 'ydir', 'normal' )
title('sum along columns projection, imag part')
xlabel( 'pages, z' )
ylabel( 'rows, y' )
daspect( [ 1, 1, 1 ] )
colorbar

subplot(233); 
imagesc( imag( squeeze( sum( expt.sample.indxrefr_data_raw, 3 ) / expt.sz.p )));
set( gca, 'ydir', 'normal' )
title('sum along pages projection, imag part')
xlabel( 'columns, x' )
ylabel( 'rows, y' )
daspect( [ 1, 1, 1 ] )
colorbar

subplot(234); 
imagesc( real( squeeze( sum( expt.sample.indxrefr_data_raw, 1 ) / expt.sz.r )));
set( gca, 'ydir', 'normal' )
title('sum along rows projection, real part')
xlabel( 'pages, z' )
ylabel( 'columns, x' )
daspect( [ 1, 1, 1 ] )
colorbar

subplot(235); 
imagesc( real( squeeze( sum( expt.sample.indxrefr_data_raw, 2 ) / expt.sz.c )));
set( gca, 'ydir', 'normal' )
title('sum along columns projection, real part')
xlabel( 'pages, z' )
ylabel( 'rows, y' )
daspect( [ 1, 1, 1 ] )
colorbar

subplot(236); 
imagesc( real( squeeze( sum( expt.sample.indxrefr_data_raw, 3 ) / expt.sz.p )));
set( gca, 'ydir', 'normal' )
title('sum along pages projection, real part')
xlabel( 'columns, x' )
ylabel( 'rows, y' )
daspect( [ 1, 1, 1 ] )
colorbar

%========

tmp1 = imag( expt.sample.indxrefr_data_raw );
tmp2 = real( expt.sample.indxrefr_data_raw );

figure; 

subplot(231); 
imagesc( squeeze( tmp1( round( 0.5 * end + 1 ), :, : )));
title('central slice, imag part')
xlabel( 'pages, z' )
ylabel( 'columns, x' )
daspect( [ 1, 1, 1 ] )
colorbar

subplot(232); 
imagesc( imag( squeeze( expt.sample.indxrefr_data_raw )));
title('sum along columns projection, imag part')
xlabel( 'pages, z' )
ylabel( 'rows, y' )
daspect( [ 1, 1, 1 ] )
colorbar

subplot(233); 
imagesc( imag( squeeze( expt.sample.indxrefr_data_raw )))
title('sum along pages projection, imag part')
xlabel( 'columns, x' )
ylabel( 'rows, y' )
daspect( [ 1, 1, 1 ] )
colorbar

subplot(234); 
imagesc( real( squeeze( expt.sample.indxrefr_data_raw )));
title('sum along rows projection, real part')
xlabel( 'pages, z' )
ylabel( 'columns, x' )
daspect( [ 1, 1, 1 ] )
colorbar

subplot(235); 
imagesc( real( squeeze( expt.sample.indxrefr_data_raw )));
title('sum along columns projection, real part')
xlabel( 'pages, z' )
ylabel( 'rows, y' )
daspect( [ 1, 1, 1 ] )
colorbar

subplot(236); 
imagesc( real( squeeze( expt.sample.indxrefr_data_raw )))
title('sum along pages projection, real part')
xlabel( 'columns, x' )
ylabel( 'rows, y' )
daspect( [ 1, 1, 1 ] )
colorbar



%=======================================================
% convert from refractive index to transmission function
%=======================================================

expt.energy = energy;                           % x-ray energy ( keV )    
expt.lambda = ( 12.4 / expt.energy ) * 1e-10;   % wavelength ( in meters )


expt.sample.TFvol = exp( +1 * ( 2 * pi / expt.lambda ) * expt.csys.z2.dLz * imag( expt.sample.indxrefr_data_raw ) ) .* ...
                    exp( -1i * ( 2 * pi / expt.lambda ) * expt.csys.z2.dLz * real( expt.sample.indxrefr_data_raw ) );
    

figure; 
subplot(131); 
imagesc( abs( squeeze( sum( expt.sample.TFvol, 1 ) / expt.sz.r )));
title('sum along rows projection')
xlabel( 'pages, z' )
ylabel( 'columns, x' )
daspect( [ 1, 1, 1 ] )

subplot(132); 
imagesc( abs( squeeze( sum( expt.sample.TFvol, 2 ) / expt.sz.c )));
title('sum along columns projection')
xlabel( 'pages, z' )
ylabel( 'rows, y' )
daspect( [ 1, 1, 1 ] )

subplot(133); 
imagesc( abs( squeeze( sum( expt.sample.TFvol, 3 ) / expt.sz.p )))
title('sum along pages projection')
xlabel( 'columns, x' )
ylabel( 'rows, y' )
daspect( [ 1, 1, 1 ] )

%====================================================================================================================================================





%=============
% create probe
%=============

expt.probe.gaussian_probe_sigma = single( probe_sigma );

expt.probe.phi = make_2Dgaussian( expt.sz.sz, 0.5 * expt.sz.sz + 1, expt.probe.gaussian_probe_sigma );

%================================
% coordinate system length scales
%================================

expt.csys.z2.Ly = expt.csys.z2.dLy * expt.sz.sz( 1 );
expt.csys.z2.Lx = expt.csys.z2.dLx * expt.sz.sz( 2 );
expt.csys.z2.Lz = expt.csys.z2.dLz * expt.sz.sz( 3 );



tic;
A = fft2( expt.probe.phi );
t1 = toc;

B = padarray( expt.probe.phi, [ 1, 1 ], 'post' );

tic;
C = fft2( B );
t2 = toc;

t1 / t2

























