



close all; clear;

wx = 15;
wy = 13;



rho = zeros( wy, wx );

% rho = rho( : );
% numelements = round( 0.1 * wx * wy );
% indices = randperm( length( rho ));
% indices = indices( 1 : numelements );
% rho( indices ) = ( 2 * rand( numelements, 1 ) - 1 ) + 1i * ( 2 * rand( numelements, 1 ) - 1 );
% f = reshape( rho, [ wy, wx ] );

rho = (1 : wx * wy);
rho = reshape( rho, [wx, wy ] ).'
% rho = reshape( rho, [wy, wx] );

f = rho;



% img_dir = '~/Documents/MATLAB/Code/data/simulated/sample/';
% f = imread([ img_dir, 'siemens-spoke.ppm' ] );
% % f = imread([ img_dir, 'airforceTP.ppm' ] );
% 
% f = rgb2hsv( f );
% f = f( :, :, 3 ) .* exp( 1i * 2 * pi * 0.0 * f( :, :, 1 ));      % brightness --> abs, hue --> phase
% f = single( imresize( f, [ 30, 20 ] ));
% 
% f = modulus_limits_scale( f, [ 0, 1 ] );
% % f = phase_limits_scale( f, 2 * pi * [ -2, 2 ] );
% 
% rho = f;
% wx = size( f, 2 );
% wy = size( f, 1 );
% 




figure; 
imagescHSV( rho ); 
daspect([1 1 1])


rho = rho( : );

Dax = zeros( wx*wy, wx*wy, 'single' );
Dax = spalloc( wx*wy, wx*wy, 2 * wx*wy );

n = wx * wy;
e = ones( n, 1 );

Dax = spdiags( e, 0, n, n );
% spy(Dax)
Dax = Dax + spdiags( -e, +wy - 0, n, n );
% spy(Dax)
Dax = Dax + spdiags( -e, wy - wy * wx - 0, n, n );
% spy(Dax)


Day = spdiags( e, 0, wy, wy );
figure; spy(Day)
Day = Day + spdiags( -e, +1, wy, wy );
figure; spy(Day)
Day = Day + spdiags( -e, 1 - wy, wy, wy );
figure; spy(Day)


ACell = repmat({Day}, 1, wx);
Day = blkdiag(ACell{:});





figure; imagesc( Dax )
figure; imagesc( Day )
figure; imagesc( rho )

Ca = Dax * double( rho );
Ca = reshape( Ca, [ wy, wx ] );

figure; 
imagesc( f )
daspect([1 1 1])

figure; 
imagesc( Ca )
daspect([1 1 1])





Ca = Day * double( rho );
Ca = reshape( Ca, [ wy, wx ] );

figure; 
imagesc( f )
daspect([1 1 1])

figure; 
imagesc( Ca )
daspect([1 1 1])













v = ones( wx*wy, 1 );
Dax = diag( v, 0 );


n = 9;
e = ones(n,1);
A = spdiags( [e -2*e e], -1:1, n, n );
full(A)


clear
wx = 5;
wy = 6;
n = wx * wy;
e = ones( n, 1 );
A = spdiags( e, 0, n, n );
full(A)
A = A + spdiags( -e, +wx, n, n );
full(A)
A = A + spdiags( -e, wx - wy * wx, n, n );
full(A)















A = sprand(1000,1000,0.005);
B = sprand(1000,1000,0.005);
Af = full(A);
Bf = full(B);

timeit(@() Af*Bf)
timeit(@() A*B)