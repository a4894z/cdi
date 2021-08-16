function [ rhof ] = fresnelpropspectrum( rhoi, p )

if p.zif == 0, rhof = rhoi; return; end
    
%==================================================================================================

% 2D coordinate system indices:

sz = single( size( rhoi ));

N = sz( 1 ) * sz( 2 );

x2 = ( -0.5 * sz( 2 ) : 0.5 * sz( 2 ) - 1 ).^2; 
y2 = transpose( ( -0.5 * sz( 1 ) : 0.5 * sz( 1 ) - 1 )).^2; 

%==================================================================================================

% x2 = gpuArray( x2 ); 
% y2 = gpuArray( y2 ); 

%==================================================================================================

% propagate from initial plane zi to final plane zf using spectrum propagation method:

temp1 = -1i * pi * p.lambda * p.zif;

temp1 = fft2( fftshift( rhoi )) .* fftshift( exp( temp1 * y2 / p.zi.Lr ^ 2 ) * exp( temp1 * x2 / p.zi.Lc ^ 2 )) / N; 
rhof = fftshift( ifft2( temp1 )) * N;
% rhof = 1i * fftshift( ifft2( temp1 ));

%==================================================================================================















% [ x, y ] = meshgrid( -0.5 * sz(2) : 0.5 * sz(2) - 1, -0.5 * sz(1) : 0.5 * sz(1) - 1 );
% 
% temp1 = -1i * pi * lambda * z;
% temp1 = fft2( fftshift( rhoi )) .* fftshift( exp( temp1 * y.^2 / Li.r ^ 2 ) .* exp( temp1 * x.^2 / Li.c ^ 2 )); 
% rhof2 = 1i * fftshift( ifft2( temp1 ));

