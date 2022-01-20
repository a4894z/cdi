%
%{

%========

restoredefaultpath; 

% addpath( genpath( '~/Documents/Science/Matlab/Code/cdi/' ));
addpath( genpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi' ));
addpath( genpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/tomolamino/testing' ));


addpath( genpath( pwd ));

%========

clear; close all; testing_multislice_actual_vs_virtual_probe

%========

%}

%=======================================================
% load sample volume, etc for multislice forward problem
%=======================================================

% load /net/s8iddata/export/8-id-ECA/Analysis/atripath/Petra_aprime_0p4sampit.mat
load /net/s8iddata/export/8-id-ECA/Analysis/atripath/ash_multislice_peco_data/Petra_aprime_0p4sampit.mat

%================================
% coordinate system length scales
%================================

resolution = resolution * 1e-6;   % um to m

expt.csys.z2.dLy = single( resolution( 2 ));
expt.csys.z2.dLx = single( resolution( 3 ));
expt.csys.z2.dLz = single( resolution( 1 ));

%==================
% energy/wavelength
%==================

expt.energy = energy;                           % x-ray energy ( keV )    
expt.lambda = ( 12.4 / expt.energy ) * 1e-10;   % wavelength ( in meters )

%========

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

expt.sample.indxrefr_data_raw = permute( data, [ 3, 2, 1 ]);
expt.sample.indxrefr_data_raw = permute( expt.sample.indxrefr_data_raw, [ 2, 1, 3 ]);

expt.sample.indxrefr_data_raw = flip( expt.sample.indxrefr_data_raw, 1 );  % = free space on top, substrate on bottom
            
%==========================================
% odd length arrays suck, modify to be even
%==========================================

% expt.sample.indxrefr_data_raw = padarray( expt.sample.indxrefr_data_raw, [ 1, 1, 0 ], 0, 'post' );

expt.sample.indxrefr_data_raw( end + 1, :, : ) = expt.sample.indxrefr_data_raw( end, :, : );
expt.sample.indxrefr_data_raw( :, end + 1, : ) = expt.sample.indxrefr_data_raw( :, end, : );

%========
            
% V = imresize3( expt.sample.indxrefr_data_raw, 0.25 );
% V = V - min( V( : ));
% 
% figure; 
% params.isoslvl = 1e-6;
% params.alphalvl = 0.35;
% isosurface_phase_surface_nopermute( abs( V ), params ); % colormap gray
% daspect([1 1 1])
% axis tight

%=================
% get problem size
%=================

expt.sz.sz      = single( size( expt.sample.indxrefr_data_raw ));
expt.sz.r       = expt.sz.sz( 1 );
expt.sz.c       = expt.sz.sz( 2 );
expt.sz.p       = expt.sz.sz( 3 );
expt.sz.sqrt_rc = sqrt( expt.sz.r * expt.sz.c );  

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

subplot(121); 
imagesc( squeeze( tmp1( :, :, round( 0.5 * end + 1 ) )));
set( gca, 'ydir', 'normal' )
title('central slice, imag part')
xlabel( 'columns, x' )
ylabel( 'rows, y' )
daspect( [ 1, 1, 1 ] )
colorbar

subplot(122); 
imagesc( squeeze( tmp2( :, :, round( 0.5 * end + 1 ) )));
set( gca, 'ydir', 'normal' )
title('central slice, real part')
xlabel( 'columns, x' )
ylabel( 'rows, y' )
daspect( [ 1, 1, 1 ] )
colorbar

%=======================================================
% convert from refractive index to transmission function
%=======================================================

% expt.sample.TFvol = exp( +1 *  ( 2 * pi / expt.lambda ) * expt.csys.z2.dLz * imag( expt.sample.indxrefr_data_raw ) ) .* ...
%                     exp( -1i * ( 2 * pi / expt.lambda ) * expt.csys.z2.dLz * real( expt.sample.indxrefr_data_raw ) );
                
                
expt.sample.TFvol = exp( -1i * expt.sample.indxrefr_data_raw * ( 2 * pi / expt.lambda ) * expt.csys.z2.dLz );
    
%========

tmp1 = abs( expt.sample.TFvol );
tmp2 = angle( expt.sample.TFvol );

figure; 

subplot(121); 
imagesc( squeeze( tmp1( :, :, round( 0.5 * end + 1 ) )));
set( gca, 'ydir', 'normal' )
title('central slice, abs')
xlabel( 'columns, x' )
ylabel( 'rows, y' )
daspect( [ 1, 1, 1 ] )
colorbar

subplot(122); 
imagesc( squeeze( tmp2( :, :, round( 0.5 * end + 1 ) )));
set( gca, 'ydir', 'normal' )
title('central slice, phase')
xlabel( 'columns, x' )
ylabel( 'rows, y' )
daspect( [ 1, 1, 1 ] )
colorbar

%========

% figure; 
% subplot(131); 
% imagesc( abs( squeeze( sum( expt.sample.TFvol, 1 ) / expt.sz.r )));
% set( gca, 'ydir', 'normal' )
% title('sum along rows projection')
% xlabel( 'pages, z' )
% ylabel( 'columns, x' )
% daspect( [ 1, 1, 1 ] )
% 
% subplot(132); 
% imagesc( abs( squeeze( sum( expt.sample.TFvol, 2 ) / expt.sz.c )));
% set( gca, 'ydir', 'normal' )
% title('sum along columns projection')
% xlabel( 'pages, z' )
% ylabel( 'rows, y' )
% daspect( [ 1, 1, 1 ] )
% 
% subplot(133); 
% imagesc( abs( squeeze( sum( expt.sample.TFvol, 3 ) / expt.sz.p )));
% set( gca, 'ydir', 'normal' )
% title('sum along pages projection')
% xlabel( 'columns, x' )
% ylabel( 'rows, y' )
% daspect( [ 1, 1, 1 ] )

%===========================================================================
% Peco's numbers for index of refraction are dodgy, guessed corrections here
%===========================================================================

expt.sample.TFvol( abs( expt.sample.TFvol ) > 1 ) = 1.0;

%=============
% create probe
%=============

expt.probe.gaussian_probe_sigma = single( [ probe_sigma( 2 ) * 1e-6 / expt.csys.z2.dLy, probe_sigma( 1 ) * 1e-6 / expt.csys.z2.dLx ] );

expt.probe.phi = make_2Dgaussian( expt.sz.sz, 0.5 * expt.sz.sz + 1, expt.probe.gaussian_probe_sigma );

expt.probe.phi = expt.probe.phi * sqrt( 1e3 / norm( expt.probe.phi, 'fro' ) ^ 2 );

%================================================
% critical distance past which we'll have alasing
%================================================

z0pad = round( 1.0 * [ expt.sz.r, expt.sz.c ] );
surrounding_zeros_r = z0pad( 1 );
surrounding_zeros_c = z0pad( 2 );

zcrit_z0zj( 1 ) = expt.csys.z2.dLy ^ 2 * surrounding_zeros_r / expt.lambda;
zcrit_z0zj( 2 ) = expt.csys.z2.dLx ^ 2 * surrounding_zeros_c / expt.lambda;

expt.csys.multslice_zcrit_z0zj = 1 * min( zcrit_z0zj );
    
%===================
% shift sample/probe
%===================

% expt.probe.phi = nocircshift2D( expt.probe.phi, [ -30, 0 ] );

%================================
% coordinate system length scales
%================================

expt.csys.z2.Ly = expt.csys.z2.dLy * expt.sz.sz( 1 );
expt.csys.z2.Lx = expt.csys.z2.dLx * expt.sz.sz( 2 );
expt.csys.z2.Lz = expt.csys.z2.dLz * expt.sz.sz( 3 );

expt.cys.z2.z_px_aliasing = floor( expt.csys.multslice_zcrit_z0zj / expt.csys.z2.dLz );

if expt.csys.z2.Lz > expt.csys.multslice_zcrit_z0zj
   
    warning( num2str( expt.cys.z2.z_px_aliasing, 'Aliasing present past z pixel = %d !)'))
    
end

%=========================
% set up multislice params
%=========================

multslice_dL = expt.csys.z2.dLz;
llambda      = expt.lambda;
Lr           = expt.csys.z2.Ly;
Lc           = expt.csys.z2.Lx;
sz           = expt.sz.sz;
% N            = expt.sz.sqrt_rc;          
         

rr2 = ( -0.5 * sz( 1 ) : 0.5 * sz( 1 ) - 1 ) .^ 2; 
cc2 = ( -0.5 * sz( 2 ) : 0.5 * sz( 2 ) - 1 ) .^ 2;   

rr2 = fftshift( transpose( rr2 ));
cc2 = fftshift( cc2 );

spectrum_prop_curv = exp( +1i * pi * llambda * multslice_dL * rr2 / Lr ^ 2 ) * exp( +1i * pi * llambda * multslice_dL * cc2 / Lc ^ 2 );


% figure; imagesc( fftshift( angle(spectrum_prop_curv) ))


% x_range = ( ( -0.5 * sz( 2 ) + 1 ) : ( 0.5 * sz( 2 ) - 0 ) );
% y_range = ( ( -0.5 * sz( 1 ) + 1 ) : ( 0.5 * sz( 1 ) - 0 ) );
% [ cc, rr ] = meshgrid( x_range, y_range );
% 
% dqx = 1 / Lc;
% dqy = 1 / Lr;
% 
% ux2_plus_uy2 = ( cc * dqx ) .^ 2 + ( rr * dqy ) .^ 2;
% spectrum_prop_curv = fftshift( exp( -1i * 2 * pi * ( multslice_dL / llambda ) * sqrt( 1 - ( llambda ^ 2 ) * ux2_plus_uy2 ))); 

clear( 'cc2', 'rr2', 'Lc', 'Lr', 'llambda', 'multslice_dL' )

%=====================
% send data TO the GPU
%=====================

sol.use_gpu = logical( 1 );  %#ok<LOGL>

if sol.use_gpu == true
    
    TF3D               = gpuArray( expt.sample.TFvol );
    phi                = gpuArray( expt.probe.phi );
    spectrum_prop_curv = gpuArray( spectrum_prop_curv );

%     supports.freespace3Dt = gpuArray( supports.freespace3Dt );
%     supports.freespace3Dg = gpuArray( supports.freespace3Dg );
%     supports.substrate3D  = gpuArray( supports.substrate3D );
%     supports.substrateTF  = gpuArray( supports.substrateTF );
    
%     cor                = gpuArray( cor );
%     aor                = gpuArray( aor );
%     rot_angle          = gpuArray( rot_angle );
%     rotvol_sz          = gpuArray( rotvol_sz );
%     vframe_r           = gpuArray( vframe_r );
%     vframe_c           = gpuArray( vframe_c );

end

%==================================================================================
% perform multislice FORWARD PROBLEM using spectrum propagation to get the exitwave
%==================================================================================

close all; 

[ phi3Df, exwv ] = multislice3Dforward( TF3D, phi, spectrum_prop_curv, sol.use_gpu );














% tic;
% A = fft2( expt.probe.phi );
% t1 = toc;
% 
% B = padarray( expt.probe.phi, [ 1, 1 ], 'post' );
% 
% tic;
% C = fft2( B );
% t2 = toc;
% 
% t1 / t2
% 
























