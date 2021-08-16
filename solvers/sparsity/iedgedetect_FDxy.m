function [ V ] = iedgedetect_FDxy( V_FDxy, FDxy )

%================================================

V_fft = ( FDxy.Dx_fft_conj .* fft2( V_FDxy.x ) + FDxy.Dy_fft_conj .* fft2( V_FDxy.y ) ) / FDxy.sqrt_rc;
% V_fft = ( FDxy.Dx_fft_conj .* fft2( V_FDxy.x ) + FDxy.Dy_fft_conj .* fft2( V_FDxy.y ));

%================================================

V_fft = V_fft .* FDxy.fft_scaling;

%================================================

if ~isfield( V_FDxy, 'q0' ) 
    
    V_fft( 1, 1 ) = 0;
    
else
    
    V_fft( 1 ,1 ) = V_FDxy.q0; 
    
end

%================================================

V = ifft2( V_fft ) * FDxy.sqrt_rc;
% V = ifft2( V_fft );

%================================================

end