function [ hpf_me, high_pass ] = hpf_blurcirc( hpf_me, diameter_pixels_x, diameter_pixels_y, ssigma_x, ssigma_y, N )

%
% High pass filter using blurred circle filter in Fourier space
%

high_pass = make_blurred_ellipsoid( diameter_pixels_x, diameter_pixels_y, ssigma_x, ssigma_y, N );
high_pass = high_pass / max(max( high_pass ));

high_pass = 1 - fftshift(high_pass);

%high_pass(high_pass == 0) = 1e-3;

%{
high_pass_zero = (high_pass == 0);

high_pass(high_pass_zero) = inf;

min_high_pass = min(high_pass(:));

high_pass(high_pass_zero) = min_high_pass * 0.1;
%}

%hpf_me = fftshift( ifft2( fft2( fftshift(hpf_me) ) .* high_pass ) );
hpf_me = ifft2( fft2( hpf_me ) .* high_pass );
