%
%{

clear; close all; testing_3dsample_from_2dproj3dsupp

%}

%====================================================================================================================================================

code_path = '~/Documents/Science/Matlab/Code/cdi';
cd( code_path )

restoredefaultpath; 
addpath( genpath( '~/Documents/Science/Matlab/Code/misc/phantom3d' ));
addpath( genpath( '~/Documents/Science/Matlab/Code/tomolamino/bresenham_line3d' ));
addpath( genpath( '~/Documents/Science/Matlab/Code/tomolamino/testing' ));

addpath( genpath( code_path ));

clear RESTOREDEFAULTPATH_EXECUTED code_path

%====================================================================================================================================================

Nr = 512;
Nc = 64;
Np = 512;

sz = [ Nr, Nc, Np ];
N = prod( sz );

aor = single( [ 0, 1, 0 ]);
aor = aor / norm( aor, 2 );

cor = ( 0.5 * sz + 1 ) + [ 0, 0, 0 ];

% 'linear' (default) | 'nearest' | 'cubic' | 'spline' | 'makima'

% interp_type = 'nearest';
interp_type = 'linear';
% interp_type = 'cubic';

%====================================================================================================================================================
% sample that lives on 2d plane slicing 3d space

V2 = single( phantom( 128 ));
V2 = imresize( V2, [ Nr, Nc ] );

abs_V_vec = abs( V2( : ));
V2 = V2 - min( abs_V_vec );
V2 = V2 / max( abs_V_vec );
V2 = abs( V2 );

V2 = exp( -1 * 0.3 * V2 ) .* exp( 2i * pi * 0.3 * V2 );

%========

V3 = ones( Nr, Nc, Np, 'single' );

% z_slice = +1 + round( 0.5 * Np + 1 );
% V3( :, :, z_slice ) = V2;

z_slice = 0 + round( 0.5 * Np + 1 );
V3( :, :, z_slice ) = V2;

% z_slice = -1 + round( 0.5 * Np + 1 );
% V3( :, :, z_slice ) = V2;

V3i = V3;

%========

S3 =  0 * V3;

% z_slice = +1 + round( 0.5 * Np + 1 );
% S3( :, :, z_slice ) = true;

z_slice = 0 + round( 0.5 * Np + 1 );
S3( :, :, z_slice ) = true;

% z_slice = -1 + round( 0.5 * Np + 1 );
% S3( :, :, z_slice ) = true;

%========

% figure;
% montage( abs( V3 ))
% colormap turbo
% drawnow
% 
% volplot.sample.isoslvl  = 0.005;
% volplot.sample.alphalvl = 0.5;
% volplot.sample.view_azel = [ 50, 30 ];
% volplot.sample.cor = 0.5 * sz + 1;
% volplot.sample.aor = [ 0, 1, 0 ];
% 
% figure; 
% set( gcf, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] );
% tomo_multislice_plotsample( (1 - abs( V3 )) .* exp(1i * angle( V3 )), volplot.sample );
% % figure; tomo_projections_along_xyz( expt.sample.vol, expt.cm.blue_light_jet );
% % view( 50, 30 )
% % % view(90, 0)

%========

% rotvol.sz  = sz;
% rotvol.cor = cor;
% rotvol.aor = aor;
% 
% [ V, ind, Xr, Yr, Zr ] = rotate_translate_volume_indices( V3, rotvol, [ 0, 0, 0, pi / 4 ] );
% 
% figure; 
% set( gcf, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] );
% tomo_multislice_plotsample( (1 - abs( V )) .* exp(1i * angle( V )), volplot.sample );
% % figure; tomo_projections_along_xyz( expt.sample.vol, expt.cm.blue_light_jet );
% % view( 50, 30 )
% % % view(90, 0)

%====================================================================================================================================================

[ C3, R3, P3 ] = meshgrid( 1 : Nc, 1 : Nr, 1 : Np );
[ C3, R3, P3 ] = deal( C3(:), R3(:), P3(:) );

RCP = transpose( [ R3, C3, P3 ]);

RCP1 = [ RCP; ones( 1, N, 'single' ) ];

T = create_3Drotmatrix( 0.9 * pi / 2, cor, aor );

RCP1rot =  T * RCP1;

Xr = RCP1rot( 2, : );
Yr = RCP1rot( 1, : );
Zr = RCP1rot( 3, : );  

Xr = reshape( Xr, sz );
Yr = reshape( Yr, sz );
Zr = reshape( Zr, sz );

%========



V3 = interp3( V3, Xr, Yr, Zr, interp_type, 1 );

% V3( isnan( V3 )) = 1;

S3 = interp3( S3, Xr, Yr, Zr, interp_type, 0 );

%========

V3( abs( V3 ) > 1 ) = 1;
S3 = round( S3 );

% abs_V_vec = abs( V3( : ));
% V3 = V3 - min( abs_V_vec );
% V3 = V3 / max( abs_V_vec );
% 
% abs_S_vec = abs( S3( : ));
% S3 = S3 - min( abs_S_vec );
% S3 = S3 / max( abs_S_vec );

%========

% figure;
% montage( abs( V3 ))
% colormap turbo
% drawnow

volplot.sample.isoslvl  = 0.005;
volplot.sample.alphalvl = 0.5;
volplot.sample.view_azel = [ 50, 30 ];
volplot.sample.cor = 0.5 * sz + 1;
volplot.sample.aor = [ 0, 1, 0 ];

figure; 
set( gcf, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] );
tomo_multislice_plotsample( (1 - abs( V3i )) .* exp( 1i * 1 * angle( V3i )), volplot.sample );
% figure; tomo_projections_along_xyz( expt.sample.vol, expt.cm.blue_light_jet );
% view( 50, 30 )
% % view(90, 0)

figure; 
set( gcf, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] );
tomo_multislice_plotsample( (1 - abs( V3 )) .* exp( 1i * 1 * angle( V3 )), volplot.sample );
% figure; tomo_projections_along_xyz( expt.sample.vol, expt.cm.blue_light_jet );
% view( 50, 30 )
% % view(90, 0)

figure; 
set( gcf, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] );
tomo_multislice_plotsample( S3, volplot.sample );
% figure; tomo_projections_along_xyz( expt.sample.vol, expt.cm.blue_light_jet );
% view( 50, 30 )
% % view(90, 0)

%====================================================================================================================================================
% true 3d sample

% V3 = single( phantom3d( 64 ));
% 
% V3 = imresize3( V3, [ 64, 64, 128 ] );
% 
% abs_V_vec = abs( V3( : ));
% V3 = V3 - min( abs_V_vec );
% V3 = V3 / max( abs_V_vec );
% 
% V3 = exp( -1 * 0.1 * V3 ) .* exp( 2i * pi * 0.1 * V3 );
% 
% figure;
% montage( abs( V3 ))
% colormap turbo
% drawnow

%====================================================================================================================================================

V2sum  = sum( V3, 3 ) / size( V3, 3 );
V2prod = prod( V3, 3 );

figure;
subplot(131);
imagescHSV( V2sum )
daspect([1 1 1])
c1 = subplot(132);
imagesc( abs( V2sum ))
daspect([1 1 1])
colorbar
set( c1, 'colormap', turbo )
c2 = subplot(133);
imagesc( angle( V2sum ))
daspect([1 1 1])
colorbar
set( c2, 'colormap', hsv )

figure;
subplot(131);
imagescHSV( V2prod )
daspect([1 1 1])
c1 = subplot(132);
imagesc( abs( V2prod ))
daspect([1 1 1])
colorbar
set( c1, 'colormap', turbo )
c2 = subplot(133);
imagesc( angle( V2prod ))
daspect([1 1 1])
colorbar
set( c2, 'colormap', hsv )

%====================================================================================================================================================

V3prod_repmat = repmat( V2prod, 1, 1, Np );
V3sum_repmat  = repmat( V2sum, 1, 1, Np );

volplot.sample.isoslvl  = 0.05;
volplot.sample.alphalvl = 0.5;
volplot.sample.view_azel = [ 50, 30 ];
volplot.sample.cor = 0.5 * sz + 1;
volplot.sample.aor = [ 0, 1, 0 ];

figure; 
set( gcf, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] );
tomo_multislice_plotsample( (1 - abs( V3prod_repmat )), volplot.sample );
% figure; tomo_projections_along_xyz( expt.sample.vol, expt.cm.blue_light_jet );
% view( 50, 30 )
% % view(90, 0)


volplot.sample.isoslvl  = 0.001;
volplot.sample.alphalvl = 0.5;
volplot.sample.view_azel = [ 50, 30 ];
volplot.sample.cor = 0.5 * sz + 1;
volplot.sample.aor = [ 0, 1, 0 ];

figure; 
set( gcf, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] );
tomo_multislice_plotsample( (1 - abs( V3sum_repmat )), volplot.sample );
% figure; tomo_projections_along_xyz( expt.sample.vol, expt.cm.blue_light_jet );
% view( 50, 30 )
% % view(90, 0)


%====================================================================================================================================================


V3prod_repmat_supp = V3prod_repmat .* S3;
V3sum_repmat_supp  = V3sum_repmat  .* S3;

V3prod_repmat_supp( V3prod_repmat_supp == 0 ) = 1;
V3sum_repmat_supp( V3sum_repmat_supp == 0 )   = 1;



volplot.sample.isoslvl  = 0.05;
volplot.sample.alphalvl = 0.5;
volplot.sample.view_azel = [ 50, 30 ];
volplot.sample.cor = 0.5 * sz + 1;
volplot.sample.aor = [ 0, 1, 0 ];

figure; 
set( gcf, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] );
tomo_multislice_plotsample( (1 - abs( V3prod_repmat_supp )), volplot.sample );
% figure; tomo_projections_along_xyz( expt.sample.vol, expt.cm.blue_light_jet );
% view( 50, 30 )
% % view(90, 0)

volplot.sample.isoslvl  = 0.001;
volplot.sample.alphalvl = 0.5;
volplot.sample.view_azel = [ 50, 30 ];
volplot.sample.cor = 0.5 * sz + 1;
volplot.sample.aor = [ 0, 1, 0 ];

figure; 
set( gcf, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] );
tomo_multislice_plotsample( (1 - abs( V3sum_repmat_supp )), volplot.sample );
% figure; tomo_projections_along_xyz( expt.sample.vol, expt.cm.blue_light_jet );
% view( 50, 30 )
% % view(90, 0)

%====================================================================================================================================================

T = create_3Drotmatrix( -0.90 * pi / 2, cor, aor );

RCP1rot =  T * RCP1;

Xr = RCP1rot( 2, : );
Yr = RCP1rot( 1, : );
Zr = RCP1rot( 3, : );  

Xr = reshape( Xr, sz );
Yr = reshape( Yr, sz );
Zr = reshape( Zr, sz );

V3f = interp3( V3prod_repmat_supp, Xr, Yr, Zr, interp_type, 1 );

volplot.sample.isoslvl  = 0.05;
volplot.sample.alphalvl = 0.5;
volplot.sample.view_azel = [ 50, 30 ];
volplot.sample.cor = 0.5 * sz + 1;
volplot.sample.aor = [ 0, 1, 0 ];

figure; 
set( gcf, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] );
tomo_multislice_plotsample( (1 - abs( V3f )) .* exp(1i * angle( V3f )), volplot.sample );
% figure; tomo_projections_along_xyz( expt.sample.vol, expt.cm.blue_light_jet );
% view( 50, 30 )
% % view(90, 0)

figure; 
set( gcf, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] );
tomo_multislice_plotsample( (1 - abs( V3i )) .* exp(1i * angle( V3i )), volplot.sample );
% figure; tomo_projections_along_xyz( expt.sample.vol, expt.cm.blue_light_jet );
% view( 50, 30 )
% % view(90, 0)

%====================================================================================================================================================