function [ pi_m ] = proj_modulus( exitwave, dffrctn_meas, missng_px_mask, N )

%propagate the exit wave to the detector plane
fft_in = fft2(fftshift( exitwave )) / N.sqrt_NxNy; 

%replace the modulus of this with the measurement but keep the phase
%temp1 = dffrctn_meas .* (fft_in ./ ( eps + abs(fft_in) )); 
% temp1 = dffrctn_meas .* exp( 1i * angle( fft_in )); 
temp1 = dffrctn_meas .* exp( 1i * angle( fft_in )) .* not( missng_px_mask ) + fft_in .* missng_px_mask; 




%back propagate back to sample plane
pi_m = fftshift(ifft2( temp1 )) * N.sqrt_NxNy;
