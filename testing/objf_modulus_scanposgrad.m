function [ f, grad ] = objf_modulus_scanposgrad( crnt_pos, sol, dffrctn_meas, N, alg, mssngdata, which )

%--------------------------------------------------------------------------------------------------

f = [];
grad = [];

%--------------------------------------------------------------------------------------------------

% temp0 = circshift( sol.smpl_trans, round(crnt_pos) );
% exitwv = temp0( N.vsy, N.vsx ) .* sol.probe;



temp0 = circshift( sol.smpl_trans, crnt_pos );
% temp0 = subpixelshift2D( sol.smpl_trans, crnt_pos );
exitwv = temp0( N.vsy, N.vsx ) .* sol.probe;

%--------------------------------------------------------------------------------------------------

if strcmp(which,'both') || strcmp(which,'val'), 
  
%   missng_px_mask = mssngdata.rect;
  [ f , ~ ] = objf_modulus( exitwv, dffrctn_meas, mssngdata, N, alg, 'val' );

end

%--------------------------------------------------------------------------------------------------

if strcmp(which,'both') || strcmp(which,'grad'), 
  
  [qc,qr] = meshgrid( (-0.5*N.Nxo+0) : (0.5*N.Nxo-1),(-0.5*N.Nyo+0) : (0.5*N.Nyo-1));
%   [qc,qr] = meshgrid( (-0.5*N.Nxo+1) : (0.5*N.Nxo-0),(-0.5*N.Nyo+1) : (0.5*N.Nyo-0));
%   [qc,qr] = meshgrid(1:N.Nxo,1:N.Nyo);
  
  qc = fftshift(qc);
  qr = fftshift(qr);

  %----------------------------------------------

  %propagate to detector:
  F_n = (fft2(fftshift( exitwv )) / N.sqrt_NxNy) .* not(mssngdata);
%   F_n = fft2( fftshift( exitwv )) / N.sqrt_NxNy;

  %----------------------------------------------

  %compute terms from sicairos/feinup 2008 paper:

  temp2 = fft2(fftshift( temp0 ));
  idft_term_rn = fftshift(ifft2( qr .* temp2 ));
  idft_term_rn = idft_term_rn( N.vsy, N.vsx );
  
  idft_term_cn = fftshift(ifft2( qc .* temp2 ));
  idft_term_cn = idft_term_cn( N.vsy, N.vsx );

  term1 = 0.5 * ( 1 - (dffrctn_meas ./ (abs(F_n) + 1E-9)) );
  
  term2 = imag( conj(F_n) .* fft2(fftshift( sol.probe .* idft_term_cn )) ) /  N.sqrt_NxNy;
  grad_xn = sum(sum( term1 .* term2 ));

  term2 = imag( conj(F_n) .* fft2(fftshift( sol.probe .* idft_term_rn )) ) /  N.sqrt_NxNy;
  grad_yn = sum(sum( term1 .* term2 ));

  %----------------------------------------------

  grad = [grad_yn grad_xn];
  
  grad = grad / norm( grad );

  %----------------------------------------------

end

%--------------------------------------------------------------------------------------------------
