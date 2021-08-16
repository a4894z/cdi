function [ sample ] = make_2Dsample( expt )

%==================================================================================================

% create a 2D sample transmission function which we will use with 2D probe modes
% to propagate exit waves ( using the projection approximation ) to the detector

%==================================================================================================
% image we use for the sample transfer function:
%==================================================================================================

% sample.img = [ expt.paths.rimgdata, '/sample/anl.ppm' ];
% sample.img = [ expt.paths.rimgdata, '/sample/brain2.ppm' ];
% sample.img = [ expt.paths.rimgdata, '/sample/simple_shapes_3.ppm' ];
sample.img = [ expt.paths.rimgdata, '/sample/cells4.ppm' ];
% sample.img = [ expt.paths.rimgdata, '/sample/siemens-spoke.ppm' ];
% sample.img = [ expt.paths.rimgdata, '/sample/ch_b.ppm' ];
% sample.img = [ expt.paths.rimgdata, '/sample/btvn.ppm' ];
% sample.img = [ expt.paths.rimgdata, '/sample/prozorov.ppm' ];
% sample.img = [ expt.paths.rimgdata, '/sample/cubes.ppm' ];
% sample.img = [ expt.paths.rimgdata, '/sample/xy3.ppm' ];
% sample.img = [ expt.paths.rimgdata, '/sample/airforceTP.ppm' ];
% sample.img = [ expt.paths.rimgdata, '/sample/tulips.ppm' ];

%================================================
% array size for sample:
%================================================

% sample.sz.r = single( round( expt.sz.r * 2.70 ));
% sample.sz.c = single( round( expt.sz.c * 2.70 ));
% sample.sz.r = single( 128 );
% sample.sz.c = single( 128 );
% sample.sz.r = single( 2048 );
% sample.sz.c = single( 2048 );

% sol.sample.sz.r = round( 1.01 * round( max( expt.spos.rs( :, 1 )) - min( expt.spos.rs( :, 1 )) + expt.sz.r ));
% sol.sample.sz.c = round( 1.01 * round( max( expt.spos.rs( :, 2 )) - min( expt.spos.rs( :, 2 )) + expt.sz.c ));
sample.sz.r = 1280;
sample.sz.c = 1280;

% sample.sz.r = expt.sz.r + 512 - 64 - 32;
% sample.sz.c = expt.sz.c + 512 - 64 - 32;
sample.sz.rc = sample.sz.r * sample.sz.c;
sample.sz.sqrt_rc = sqrt( sample.sz.rc );
sample.sz.sz = [ sample.sz.r, sample.sz.c ]; 

%================================================
% define params for sample transmission:
%================================================

sample.expofimgabs = false;
sample.abs_minmax = single( [ +0.6000, +0.999 ] );
sample.phs_minmax = single( 0.99 * [ -1.0, +1.0 ] * pi );
sample.phs_offset = single( 0.0 * pi );

%================================================
% make the sample:
%================================================

% load an image, create it to hsv, val <--> abs, hue <--> phs 
[ tmp1 ] = image2complex( sample.img, 3, 3 );

sample.T = imresize( tmp1, 1 * sample.sz.sz, 'box' );
% sample.T = imresize2D( tmp1, 1 * sample.sz.sz );

%========

% % resize to desired problem size for the sample transmission function:
% % sample.T = imresize( tmp1, 0.333 * sample.sz.sz );
% sample.T = imresize2D( tmp1, 0.333 * sample.sz.sz );
% % sample.T = padarray( sample.T, round( [ 0.5 * ( sample.sz.r - size( sample.T, 1 )), 0.5 * ( sample.sz.c - size( sample.T, 2 )) ] ));
% sample.T = zeropadarray( sample.T, 0.5 * round2even( [  ( sample.sz.r - size( sample.T, 1 )), ( sample.sz.c - size( sample.T, 2 )) ] ) );

%========

% % "stretch out" sample in vertical direction:
% sample.T = truncatearray( sample.T, [ 16, sample.sz.c ] );
% %sample.T = padarray( sample.T, [ 0, 256 ] );
% sample.T = imresize( sample.T, sample.sz.sz );
% sample.T = circshift( sample.T, [ 30, 15 ] );

%========

% make sample "smaller":
% sample.T = zeropadarray( sample.T, [ 128, 128 ] );
% sample.T = imresize2D( sample.T, 1 * sample.sz.sz );

%================================================

sample.sz.r = single( size( sample.T, 1 ));
sample.sz.c = single( size( sample.T, 2 ));
sample.sz.rc = sample.sz.r * sample.sz.c;
sample.sz.sqrt_rc = sqrt( sample.sz.rc );
sample.sz.sz = [ sample.sz.r, sample.sz.c ]; 

% view slices for when computing the exit wave views at a particular scan position
% sample.vs( 1, :  ) = round( ( 0.5 * ( sample.sz.r - expt.sz.r ) + 1 ) : ( 0.5 * ( sample.sz.r + expt.sz.r )));
% sample.vs( 2, :  ) = round( ( 0.5 * ( sample.sz.c - expt.sz.c ) + 1 ) : ( 0.5 * ( sample.sz.c + expt.sz.c )));
sample.vs.r = round( ( 0.5 * ( sample.sz.r - expt.sz.r ) + 1 ) : ( 0.5 * ( sample.sz.r + expt.sz.r )));
sample.vs.c = round( ( 0.5 * ( sample.sz.c - expt.sz.c ) + 1 ) : ( 0.5 * ( sample.sz.c + expt.sz.c )));

%================================================

% sample.T = exp( -1 * abs( sample.T )) .* exp( 1i * angle( sample.T ));
% sample.T = exp( -1 * abs( sample.T )) .* exp( 1i * abs( sample.T ));

%================================================

sample.T = modulus_limits_scale( sample.T, sample.abs_minmax );
sample.T = phase_limits_scale( sample.T, sample.phs_minmax );

% sample.T = abs( sample.T );

%================================================

% QUICK HACK: for an opaque region surrounding an isolated sample:

% sample.T = (abs( sample.T ) > 0.333) .* sample.T;

%================================================
