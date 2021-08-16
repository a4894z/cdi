function [ V_FDxy ] = edgedetect_FDxy( V, FDxy )

% RENAME TO CONVOLUTION EDGE DETECTION
%==================================================================================================

temp1 = fft2( V ) / FDxy.sqrt_rc;
% temp1 = fft2( V );

%=======================

% keep track of the q = 0 pixel of the fft2 of the current V. 
% this will be used in the inverse FDxy step to get the p
% proper scaling of the exit wave based on the q = 0 pixel from
% the measurement (assuming no beamstop)

V_FDxy.q0 = temp1( 1, 1 ); 

%=======================

V_FDxy.x = ifft2( FDxy.Dx_fft .* temp1 ) * FDxy.sqrt_rc;
V_FDxy.y = ifft2( FDxy.Dy_fft .* temp1 ) * FDxy.sqrt_rc;

% V_FDxy.x = ifft2( FDxy.Dx_fft .* temp1 ); 
% V_FDxy.y = ifft2( FDxy.Dy_fft .* temp1 ); 

%==================================================================================================

end

