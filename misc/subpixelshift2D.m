function y = subpixelshift2D( x, delta )
%
% use fourier shift theorem to do subpixel shifting
%

[ N, M ] = size( x );

X = fft2( x );

% floors take care of odd-length signals.
x_shift = exp( -2 * pi * 1i * delta(2) * [ 0 : floor( M / 2 ) - 1, floor( -M / 2 ) : -1 ]  / M );
y_shift = exp( -2 * pi * 1i * delta(1) * [ 0 : floor( N / 2 ) - 1, floor( -N / 2 ) : -1 ]' / N );

% enforce conjugate symmetry. otherwise this frequency component has no
% corresponding negative frequency to cancel out its imaginary part.
if mod( N, 2 ) == 0, y_shift( N / 2 + 1 ) = real( y_shift( N / 2 + 1 ) ); end 
if mod( M, 2 ) == 0, x_shift( M / 2 + 1 ) = real( x_shift( M / 2 + 1 ) ); end

Y = X .* ( y_shift * x_shift );

y = ifft2( Y );

% There should be no imaginary component (for real input signals) but due to numerical effects some remnants remain.
if isreal( x ), y = real( y ); end
