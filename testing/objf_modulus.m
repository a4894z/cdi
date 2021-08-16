function [ f , grad_f ] = objf_modulus( crrnt_view, dffrctn_meas, missng_px_mask, N, alg, which )


% crrnt_view is the current exit wave iterate (N.Ny x N.Nx single array)
% 
% cnstrnt is a struct containing the diffraction measurement constraint, the support, and a mask defining which elements of
% the diffraction measurement are non-zero
% 
% N is a struct with the sizes of the arrays used
% 
% alg is a struct with parameters used for defining and computing the error metric, see MGS Opt Express 2008
%   alg.modmetrc_gam =  0.5 usually (is what Marchesini RSI 2007 uses)
%   alg.modmetrc_del ~  eps (a small number, used for when the diffraction measurement = 0, as far as a I can tell, we don't 
%                       need and can remove this variable...is used in MGS Opt Express 2008)
% 
% which is a string telling us if we want the objective function value, its gradient, or both

f = [];
grad_f = [];

% if ~exist(alg,'var'), alg.modmetrc_gam = 0.5; end

if alg.modmetrc_gam == 0.5,
  
    pi_m = proj_modulus( crrnt_view, dffrctn_meas, missng_px_mask, N );
    temp1 = pi_m - crrnt_view;

  
%     fft_in = fft2(fftshift( exitwave )) / N.sqrt_NxNy; 
%     %temp1 = dffrctn_meas .* (fft_in ./ ( eps + abs(fft_in) )); 
%     % temp1 = dffrctn_meas .* exp( 1i * angle( fft_in )); 
%     temp1 = dffrctn_meas ; 
  

    if strcmp(which,'both') || strcmp(which,'val'), 
    f = sum(sum( abs( temp1 ).^2 )); end

    if strcmp(which,'both') || strcmp(which,'grad'),
    grad_f = -2*temp1; end

else
  
  tilde_psi = fft2( crrnt_view ) / N.sqrt_NxNy;

  temp1 = missng_px_mask .* ( abs( tilde_psi ).^2 + alg.modmetrc_del );
  temp2 =  ( temp1 .^ alg.modmetrc_gam - ( dffrctn_meas .^ 2 + alg.modmetrc_del ) .^ alg.modmetrc_gam );

  if strcmp(which,'both') || strcmp(which,'val'), 
    f = sum(sum( temp2 .^2 )); end

  if strcmp(which,'both') || strcmp(which,'grad'),
    A_q = alg.modmetrc_gam * temp2 .* ( eps + temp1 .^ (alg.modmetrc_gam - 1) );
    grad_f = 4 * ifft2( A_q .* tilde_psi ) * N.sqrt_NxNy ;      
  end

end

