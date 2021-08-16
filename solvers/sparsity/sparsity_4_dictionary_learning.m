
%{

% see http://www.numerical-tours.com/matlab/sparsity_4_dictionary_learning/#1




clear; close all; 
restoredefaultpath
cd( '~/Documents/MATLAB/Code/cdi/' );
addpath(genpath(pwd));












clear; close all; sparsity_4_dictionary_learning






% have one config file to set up an experiment, 
% but use options flags for num of scan pos, start pos on sample, 
% type of beam (pinhole, focusing optics)



%}

%==================================================================================================
%---------------------------------- load initial this and that ------------------------------------
%==================================================================================================

tmp0 = load( '~/Documents/MATLAB/Code/cdi/run/misctesting2DTPA/simPRexpt_2DTPA.mat' );

expt = tmp0.expt;
sol = tmp0.sol;

clearvars -except expt sol;

%==================================================================================================
%-------- define image we will learn the sparse analysis and synthesis representation for ---------
%==================================================================================================

%{

img_dir = '~/Documents/MATLAB/Code/data/simulated/sample/';

% f = imread([ img_dir, 'btvn.ppm' ] );
% f = imread([ img_dir, 'baboon.ppm' ] );
f = imread([ img_dir, 'fractal1.ppm' ] );
% f = imread([ img_dir, 'cells2.ppm' ] );
% f = imread([ img_dir, 'dislocationsA.ppm' ] );
% f = imread([ img_dir, 'dislocationsB.ppm' ] );
% f = imread([ img_dir, 'hsv.ppm' ] );
% f = imread([ img_dir, 'leaf2b.ppm' ] );
% f = imread([ img_dir, 'spiral_zp.ppm' ] );
% f = imread([ img_dir, 'SheppLogan2.ppm' ] );
% f = imread([ img_dir, 'tulips.ppm' ] );
% f = imread([ img_dir, 'xy.ppm' ] );
% f = imread([ img_dir, 'xy3.ppm' ] );
% f = imread([ img_dir, 'simple_shapes_3.ppm' ] );
% f = imread([ img_dir, 'spiral.ppm' ] );
% f = imread([ img_dir, 'spiral2.ppm' ] );
% f = imread([ img_dir, 'spiral3.ppm' ] );
% f = imread([ img_dir, 'siemens-spoke.ppm' ] );
% f = imread([ img_dir, 'airforceTP.ppm' ] );

%=======================

f = rgb2hsv( f );
f = f( :, :, 3 ) .* exp( 1i * 2 * pi * 0.99 * f( :, :, 1 ));      % brightness --> abs, hue --> phase

% temp1 = load( '/home/l318967/Documents/MATLAB/Code/AshPRCombined/combined/fcdi_simulated_experiment_data.mat' );
% f = temp1.sol.probe;
% clear('temp1')

% % phase dictionary learning:
% f = angle( f ) + pi;
% f = f / max( f(:) );

% f = abs( f );
% f = f / max( f(:) );


%=======================

% load( '/home/ash/Documents/MATLAB/Code/masochist_lizard/numerical_tours/sparsity_4_dictionary_learning/OEvoronoi_16.mat' );
% f = exitwaveOE;
% f = f / max( max( abs( f )));
% 
% NN = 256;
% f = f( round( ( size( f, 1 ) - NN ) / 2 + 1 ) : round( ( size( f, 1 ) + NN ) / 2), round( ( size( f, 2 ) - NN ) / 2 + 1 ) : round( ( size( f, 2 ) + NN ) / 2 ));
% clear( 'exitwaveOE', 'NN' )

%=======================

% f = single( imresize( f, [ 32, 32 ] ));
% f = single( imresize( f, [ 64, 64 ] ));
% f = single( imresize( f, [ 96, 96 ] ));
% f = single( imresize( f, [ 128, 128 ] ));
f = single( imresize( f, [ 128 + 64, 128 ] ));
% f = single( imresize( f, [ 256, 256 ] ));
% f = single( imresize( f, [ 384, 256 ] ));
% f = single( imresize( f, [ 384, 384 ] ));
% f = single( imresize( f, [ 512, 512 ] ));
% f = single( imresize( f, [ 640, 640 ] ));

f = modulus_limits_scale( f, [ 0, 1 ] );
f = phase_limits_scale( f, 2 * pi * [ -2, 2 ] );

%=======================

% % f = single( imresize( f, [ 4, 4 ] ));
% 
% f = single( imresize( f, [ 3, 3 ] ));
% f = f * 0;
% 
% % f( 1, 1 ) = 1 * exp( 1i * pi * 0.2 ); 
% % f( 2, 2 ) = 1 * exp( 1i * pi * 0.5 ); 
% % f( 3, 3 ) = 1 * exp( 1i * pi * 0.6 ); 
% % 
% % f( 2, 1 ) = 0.7 * exp( 1i * pi * 0.9 ); 
% % f( 3, 1 ) = 0.7 * exp( 1i * pi * 0.9 ); 
% % f( 3, 2 ) = 0.7 * exp( 1i * pi * 0.9 ); 
% % 
% % f( 1, 2 ) = 0.5 * exp( -1i * pi * 0.9 ); 
% % f( 1, 3 ) = 0.5 * exp( -1i * pi * 0.9 ); 
% % f( 2, 3 ) = 0.5 * exp( -1i * pi * 0.9 ); 
% 
% 
% 
% f( 2, 3 ) = 0.5 * exp( -1i * pi * 0.3 ); 
% f( 2, 1 ) = 0.5 * exp( -1i * pi * 0.3 ); 
% 
% f( 1, 2 ) = 0.5 * exp( -1i * pi * 0.3 ); 
% f( 3, 2 ) = 0.5 * exp( -1i * pi * 0.3 ); 
% 
% f( 1, 1 ) = 1.0 * exp( +1i * pi * 0.3 ); 
% f( 1, 3 ) = 1.0 * exp( +1i * pi * 0.3 ); 
% 
% f( 3, 1 ) = 1.0 * exp( +1i * pi * 0.5 ); 
% f( 3, 3 ) = 1.0 * exp( +1i * pi * 0.5 ); 
% 
% f( 2, 2 ) = 0.25 * exp( +1i * pi * 0.8 ); 

%=======================

% sigma = .06;    % noise "level" for additive gaussian
% f = f + sigma * randn( sol.sz.r );

%=======================

% fN = random( 'Poisson', 500 * abs( f )); 
% f = ( fN / 500 ) .* exp( 1 * 1i * angle( f )); 

%=======================

% figure; imagesc( abs(f), [ min( abs(f(:)) ), max( abs(f(:)) )] ); daspect([1 1 1]); colormap gray; title('Original Image')

%=======================

% image size 
sol.sz.r = size( f, 1 );       
sol.sz.c = size( f, 2 ); 
sol.sz.rc = sol.sz.r * sol.sz.c;
sol.sz.sqrt_rc = sqrt( sol.sz.rc );

%}

%==================================================================================================

[ sol.sDL ] = setup_synthesis_dictionarylearning( sol.sz );

[ sol.aDL ] = setup_analysis_dictionarylearning( sol.sz );

%==================================================================================================

ss = 5;
[ f, ~ ] = enforce_2DTPAsposview( expt.probe.P, expt.sample.TF, expt.sample.vs, expt.spos.rs( ss, : ), 'subpx' );

% [ Ys, sol.sDL ] = image2trainingsetpatches( f, sol.sDL, sol.sz  );
% [ Ya, sol.aDL ] = image2trainingsetpatches( f, sol.aDL, sol.sz  );
% [ fs ] = trainingsetpatches2image( Ys, sol.sDL, sol.sz );
% [ fa ] = trainingsetpatches2image( Ya, sol.aDL, sol.sz );
% 
% % figure; imagesc( abs(f), [ min( abs(f(:)) ), max( abs(f(:)) )] ); daspect([1 1 1 ]); colormap gray
% % figure; imagesc( abs(fs), [ min( abs(f(:)) ), max( abs(f(:)) )] ); daspect([1 1 1 ]); colormap gray
% % figure; imagesc( abs(fa), [ min( abs(f(:)) ), max( abs(f(:)) )] ); daspect([1 1 1 ]); colormap gray
% 
% figure; imagescHSV( f )
% figure; imagescHSV( fs )
% figure; imagescHSV( fa )

%=======================

% % display the dictionary.
% figure; 
% subaxis(1,1,1,'SpacingVert',0,'MR',0.01, 'ML',0.08,'MT',0.06, 'MB',0.08); 
% plot_dictionnary( sol.sDL, [ 12, 12 ] );
% % axis off;
% % daspect([1 1 1])
% title('Subset of Abs( Sythesis Starting Dictionary )')



%==================================================================================================
%------------------------------- sparse ANALYSIS learning ----------------------------------------
%==================================================================================================

%{

sol.aDL.ntot = 100;
sol.aDL.ncoef = 1;
sol.aDL.ndict = 1;
sol.aDL.E0 = [];
sol.aDL.t = [];

sol.aDL.quiet = 'no';

[ f1a, sol.aDL ] = proj_sparse_analysis_learning( f, sol.aDL, sol.sz );

fprintf( '\n' )

%=======================

limits = [ 0, 1 ];

figure; 
set( gcf, 'units', 'normalized', 'outerposition', [0 0 1.0 1.0] )

% subplot(221); 
subaxis(2,2,1,'SpacingVert',0,'MR',0.01, 'ML',0.01,'MT',0.05, 'MB',0.05); 
imagescHSV( f ); daspect([1 1 1]); title('Original Image'); axis off
% subplot(223); 
subaxis(2,2,3,'SpacingVert',0,'MR',0.01, 'ML',0.01,'MT',0.1, 'MB',0.05); 
imagesc( abs( f ), limits ); daspect([1 1 1]); colormap( expt.cm.blj ); axis off; %colorbar; 
title('| Original Image |')

% subplot(222);
subaxis(2,2,2,'SpacingVert',0,'MR',0.01, 'ML',0.01,'MT',0.05, 'MB',0.05); 
imagescHSV( f1a ); daspect([1 1 1]); title('Analysis Sparse Image'); axis off; 
% subplot(224); 
subaxis(2,2,4,'SpacingVert',0,'MR',0.01, 'ML',0.01,'MT',0.1, 'MB',0.05); 
imagesc( abs( f1a ), limits ); daspect([1 1 1]); colormap( expt.cm.blj ); axis off; %colorbar; 
title('| Analysis Sparse Image |')

erra = (( abs( f1a ) ~= 0 ) .* abs( f - f1a ));

figure; 
set( gcf, 'units', 'normalized', 'outerposition', [0 0 1.0 1.0] )
subplot(121); 
imagesc( erra ); daspect([1 1 1]); colormap( expt.cm.blj ); axis off; %colorbar; 
title('| analysis - orig |')
%title('| analysis | - | orig |')
%colorbar( 'peer', spea, 'Position', [0.849930843706778 0.1220703125 0.0165975103734443 0.0673828125], 'Color', [1 1 1] );
colorbar; 

subplot(122); 
set( gcf, 'units', 'normalized', 'outerposition', [0 0 1.0 1.0] )
semilogy( sol.aDL.E0, '-o', 'LineWidth', 2 );

%}

%==================================================================================================
%------------------------------- sparse SYNTHESIS learning ----------------------------------------
%==================================================================================================

sol.sDL.ntot = 100;
sol.sDL.ncoef = 1;
sol.sDL.ndict = 1;
sol.sDL.EC = [];
sol.sDL.ED = [];
sol.sDL.t = [];

sol.sDL.quiet = 'no';

[ f1s, sol.sDL ] = proj_sparse_synthesis_learning( f, sol.sDL, sol.sz );

fprintf( '\n' )

%=======================

% limits = [ 0, 1 ];
limits = [ 0, max( abs( f( : ))) ];

figure; 
%set( gcf, 'units', 'normalized', 'outerposition', [0 0 1.0 1.0] )

subplot(231); 
% subaxis(2,2,1,'SpacingVert',0,'MR',0.01, 'ML',0.01,'MT',0.05, 'MB',0.05); 
imagescHSV( f ); daspect( [ 1, 1, 1 ] ); title('Original Image'); axis off
subplot(234); 
% subaxis(2,2,3,'SpacingVert',0,'MR',0.01, 'ML',0.01,'MT',0.1, 'MB',0.05); 
imagesc( abs( f ), limits ); daspect( [ 1, 1, 1 ] ); colormap( expt.cm.blj ); axis off; %colorbar; 
title('| Original Image |')

subplot(232);
% subaxis( 2, 2, 2, 'SpacingVert', 0, 'MR', 0.01, 'ML', 0.01, 'MT', 0.05, 'MB', 0.05 ); 
imagescHSV( f1s ); daspect( [ 1, 1, 1 ] ); title( 'Synthesis Sparse Image' ); axis off; 
subplot(235); 
% subaxis(2,2,4,'SpacingVert',0,'MR',0.01, 'ML',0.01,'MT',0.1, 'MB',0.05); 
imagesc( abs( f1s ), limits ); daspect( [ 1, 1, 1 ] ); colormap( expt.cm.blj ); axis off; %colorbar; 
title('| Synthesis Sparse Image |')

%=======================

errs = (( abs( f1s ) ~= 0 ) .* abs( f - f1s ));

% figure; 
%set( gcf, 'units', 'normalized', 'outerposition', [0 0 1.0 1.0] )
subplot(233); 
imagesc( errs ); daspect( [ 1, 1, 1 ] ); colormap( expt.cm.blj ); axis off; %colorbar; 
title('| synthesis - orig |')
%title('| analysis | - | orig |')
%colorbar( 'peer', spea, 'Position', [0.849930843706778 0.1220703125 0.0165975103734443 0.0673828125], 'Color', [1 1 1] );
colorbar; 

subplot(236); 
%set( gcf, 'units', 'normalized', 'outerposition', [0 0 1.0 1.0] )
hold on
semilogy( sol.sDL.EC, '-o', 'LineWidth', 3, 'color', [0, 0, 0] );
semilogy( sol.sDL.ED, '-x', 'LineWidth', 2, 'color', [0, 1, 0] );
hold off

%=======================

tmp0 = sol.sDL.D;
% tmp0( ~any( tmp0, 2 ), : ) = [];  %rows
tmp0( :, ~any(tmp0,1) ) = [];  %columns

sDL2 = sol.sDL;
% sDL2.D( :, ~any(sDL2.D,1) ) = [];  %columns

figure; plot_dictionnary( sDL2, [ 4, 4 ] );

%=======================

figure; 
[ V, S, U ] = svd( f ); 
semilogy( diag( S ), '-', 'linewidth', 4 );
[ V, Ss, U ] = svd( f1s );
hold on; semilogy( diag( Ss ),'-o'); hold off; 
grid on
ylim( [ 10^-6, 10^4 ] )
legend('svd( orig image )', 'svd( synthesis sparse image )' );


return

%==================================================================================================
%--------------------------------------- plot the results  ----------------------------------------
%==================================================================================================

figure; 
set( gcf, 'units', 'normalized', 'outerposition', [0 0.6 1.0 0.4] )

% errs = (( abs( f1s ) ~= 0 ) .* abs( abs( f ) - abs( f1s )));
% erra = (( abs( f1a ) ~= 0 ) .* abs( abs( f ) - abs( f1a )));

errs = (( abs( f1s ) ~= 0 ) .* abs( f - f1s ));
erra = (( abs( f1a ) ~= 0 ) .* abs( f - f1a ));

% errs = (( abs( f1s ) ~= 0 ) .* abs( abs( f ) - abs( f1s ))) ./ ( 1e-7 + ( abs( f1s ) ~= 0 ) .* abs( f ));
% erra = (( abs( f1a ) ~= 0 ) .* abs( abs( f ) - abs( f1a ))) ./ ( 1e-7 + ( abs( f1a ) ~= 0 ) .* abs( f ));

% errs = (( abs( f1s ) ~= 0 ) .* abs( f - f1s )) ./ ( 1e-7 + ( abs( f1s ) ~= 0 ) .* abs( f ));
% erra = (( abs( f1a ) ~= 0 ) .* abs( f - f1a )) ./ ( 1e-7 + ( abs( f1a ) ~= 0 ) .* abs( f ));

% errs = (( abs( f1s ) ~= 0 ) .* abs( abs( f1s ))) ./ ( 1e-7 + ( abs( f1s ) ~= 0 ) .* abs( f ));
% erra = (( abs( f1a ) ~= 0 ) .* abs( abs( f1a ))) ./ ( 1e-7 + ( abs( f1a ) ~= 0 ) .* abs( f ));

spes = subplot(121);
imagesc( errs ); daspect([1 1 1]); colormap( expt.cm.blj ); axis off; 
title('| synthesis - orig |')
%title('| synthesis | - | orig |')
%colorbar( 'peer', spes, 'Position', [0.863070539419088 0.849609375 0.0145228215767633 0.064453125], 'Color', [1 1 1] );
colorbar; 

spea = subplot(122);
imagesc( erra ); daspect([1 1 1]); colormap( expt.cm.blj ); axis off; %colorbar; 
title('| analysis - orig |')
%title('| analysis | - | orig |')
%colorbar( 'peer', spea, 'Position', [0.849930843706778 0.1220703125 0.0165975103734443 0.0673828125], 'Color', [1 1 1] );
colorbar; 

%================================================

figure; 
set( gcf, 'units', 'normalized', 'outerposition', [0 0.0 1.0 0.4] )

%figure; 
subplot(132)
plot( 1 : sol.sDL.ntot, sol.sDL.t, 'LineWidth', 2 )
hold on
plot( 1 : sol.aDL.ntot, sol.aDL.t, 'LineWidth', 2 )
%axis tight;
hold off
grid on
legend('Timing, Synthesis Learning', 'Timing, Analysis Learning');

%figure; 
subplot(133)
[ V, S, U ] = svd( f ); 
semilogy( diag( S ), '-', 'linewidth', 4 );
[ V, Ss, U ] = svd( f1s );
hold on; semilogy( diag( Ss ),'-o'); hold off; 
[ V, Sa, U ] = svd( f1a ); 
hold on;  semilogy( diag( Sa ),'-x');  hold off;
grid on
ylim( [ 10^-6, 10^4 ] )
legend('svd( orig image )', 'svd( synthesis sparse image )', 'svd( analysis sparse image )');

subplot(131)
semilogy( sol.sDL.E0, 'LineWidth', 2 );
hold on;
semilogy( sol.aDL.E0, 'LineWidth', 2 );
hold off;
axis tight;
grid on
legend( '|| Y - Ds Cs ||_2^2', '|| Ca - Da Y ||_2^2' );

% subplot(131)
% semilogy( 1 : 2 * sol.sDL.ntot, sol.sDL.E0);
% hold on;
% semilogy( 1 : 2 : 2 * sol.sDL.ntot, sol.sDL.E0( 1 : 2 : 2 * sol.sDL.ntot ), '*', 'LineWidth', 2 );
% semilogy( 2 : 2 : 2 * sol.sDL.ntot, sol.sDL.E0( 2 : 2 : 2 * sol.sDL.ntot ), 'o', 'LineWidth', 2 );
% semilogy( 1 : 2 * sol.sDL.ntot, sol.aDL.E0 );
% semilogy( 1 : 2 : 2 * sol.sDL.ntot, sol.aDL.E0( 1 : 2 : 2 * sol.sDL.ntot ), '*', 'LineWidth', 2 );
% semilogy( 2 : 2 : 2 * sol.sDL.ntot, sol.aDL.E0( 2 : 2 : 2 * sol.sDL.ntot ), 'o', 'LineWidth', 2 );
% axis tight;
% hold off
% grid on
% legend('|| Y - Ds Cs ||_2^2', 'After coefficient update', 'After dictionary update', '|| Ca - Da Y ||_2^2', 'After coefficient update', 'After dictionary update');

%================================================


return

%{

sol.aDL.theta = sol.aDL.theta * 0;

aDL_C = sol.aDL.C( 1 : sol.aDL.Ni, : );
aDL_C = sol.aDL.C( sol.aDL.Ni + 1 : 2 * sol.aDL.Ni, : );
aDL_C = sol.aDL.C( 2 * sol.aDL.Ni + 1 : 3 * sol.aDL.Ni, : );

[ f2 ] = trainingsetpatches2image( aDL_C, sol.aDL, sol.sz );
%figure; imagesc(log10(1+abs(f2)))
figure; imagescHSV( f2 ); daspect([ 1, 1, 1 ])

%}

%{
close all

figure; imagescHSV( f ); daspect([1 1 1]); axis off;
figure; imagescHSV( f1a ); daspect([1 1 1]); axis off;


sly = 59 : 150;
slx = 2 : 40;
figure; imagescHSV( f( sly, slx ) ); daspect([1 1 1]); axis off;
figure; imagescHSV( f1a( sly, slx ) ); daspect([1 1 1]); axis off;

sly = 148 : 241;
slx = 131 : 250;
figure; imagescHSV( f( sly, slx ) ); daspect([1 1 1]); axis off;
figure; imagescHSV( f1a( sly, slx ) ); daspect([1 1 1]); axis off;
%}


%{


sparsity = edge_detect_xy_forward_diff( [], sol.sz );


[ prtl_rho ] = forward_edgedetect( f, sparsity );


max_f = max( prtl_rho.x( : ));
f_edge_abs_sort = sort( abs( prtl_rho.x( : )), 'descend');
sel = f_edge_abs_sort( round( numel( f_edge_abs_sort ) * 0.010 ));

prtl_rho.x = prtl_rho.x .* ( abs( prtl_rho.x ) >= sel );


max_f = max( prtl_rho.y( : ));
f_edge_abs_sort = sort( abs( prtl_rho.y( : )), 'descend');
sel = f_edge_abs_sort( round( numel( f_edge_abs_sort ) * 0.010 ));

prtl_rho.y = prtl_rho.y .* ( abs( prtl_rho.y ) >= sel );


[ rho ] = inverse_edgedetect( prtl_rho, sparsity, sol.sz );


figure; imagesc( abs( prtl_rho.x ) ); axis off; daspect([1 1 1])

figure; imagesc( abs( f ), [0, 1 ]) ; axis off; daspect([1 1 1])
figure; imagesc( abs( rho ), [0, 1]) ; axis off; daspect([1 1 1])

figure; imagescHSV( f ); axis off; daspect([1 1 1])
figure; imagescHSV( rho ); axis off; daspect([1 1 1])



Dx = 0 * f;
Dx( 1, 1 ) = +1;
Dx( 1, 2 ) = -1;

f_edge = ifft2( fft2( Dx ) .* fft2( f ));

max_f = max( f_edge( : ));
f_edge_abs_sort = sort( abs( f_edge( : )), 'descend');
sel = f_edge_abs_sort( round( numel( f_edge_abs_sort ) * 0.013 ));

f_edgeT = f_edge .* ( abs( f_edge ) >= sel );

f_edgeTi = ifft2( conj( fft2( Dx ) ) .* fft2( f_edgeT ) ./ ( 1e-7 + abs( fft2( Dx )) .^2  ));


figure; imagesc( abs( f_edge ) ) 
figure; imagesc( abs( f_edgeT ) ) 
figure; imagesc( abs( f_edgeTi ) ) 

%}


 



%{



% LOOK UP FORWARD DIFFERENCE CONVOLUTION MATRIX OPERATOR



Da = zeros(  sol.aDL.Na, sol.aDL.Ni, 'single' );
Da( 1, 1 ) = +1;
Da( 1, 2 ) = -1;


for ii = 2 : ( sol.aDL.w - 1 )
    tmp1 = circshift( Da( 1, : ), [ 0, ii - 1 ] );
    Da( ii, : ) = tmp1;
end


figure; imagesc( Da )


temp0 = Da * Ya;
temp0 = circshift( reshape( temp0, [sol.aDL.w,sol.aDL.w] ), [1, 0] );
figure; imagesc( abs( temp0 )); 

% COMPARE ABOVE TO CONVOLUTION


dc2 = zeros( sol.aDL.w, sol.aDL.w, 'single' );
dc2(1,1) = +1;
dc2(2,1) = -1;

temp1 = ifft2( ( fft2( dc2 ) ) .* ( fft2( f ) ) );

figure; imagesc( abs( temp1 ))

figure; imagesc( abs( temp1) - abs(temp0))

figure; imagesc( angle( temp1 ) - angle(temp0))


5;
%}


