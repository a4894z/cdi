function [ V ] = lpf_blurcirc( V, diam, stdev )

%
% High pass filter using blurred circle filter in Fourier space
%

sz = size( V );

low_pass = make_2Dellipsoid( sz, diam );

[ low_pass ] = lpf_gauss( low_pass, stdev );

low_pass =  fftshift( low_pass / max( abs( low_pass( : ))));

V = fftshift( ifft2( fft2( fftshift( V )) .* low_pass ));
