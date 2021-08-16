function img_sp_shft = subpixelshift_fft2_old( img, delta_r, delta_c )


[ nr, nc ] = size( img );

Nr = ifftshift( -fix(nr/2) : ceil(nr/2) - 1 );
Nc = ifftshift( -fix(nc/2) : ceil(nc/2) - 1 );

[Nc,Nr] = meshgrid( Nc, Nr );

img_sp_shft = ifft2( fft2(img) .* exp( 1i * 2 * pi * ( delta_r * Nr / nr + delta_c * Nc / nc ) ) ) ;

% There should be no imaginary component (for real input signals) but due to numerical effects some remnants remain.
if isreal( img ), img_sp_shft = real( img_sp_shft ); end

