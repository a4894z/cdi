function [ prtf ] = computePRTFvsQ( meas_modulus, soln_modulus, N, mssgdta )

%#ok<*BDSCI>
%-------------------------------------------------------------------------

%keep this for the time being, compare speeds for defining it here vs passing it
[xx,yy] = meshgrid(1:N.Nc,1:N.Nr);

%-------------------------------------------------------------------------

if ~exist('mssgdta','var')
  mssgdta.rect = zeros(N.Nr,N.Nc);
end

beam_stop = not(fftshift(mssgdta.rect));

%-------------------------------------------------------------------------

%the q = 0 value:
jj = 1;
o_diam_circ = (((N.Nc/2+1-xx)./jj).^2 + ((N.Nr/2+1-yy)./jj).^2 < 1);

kk = 1;
prtf(kk) = sum(sum(o_diam_circ .* soln_modulus)) / ( eps + sum(sum(o_diam_circ .* meas_modulus)) );

kk = kk + 1;

%-------------------------------------------------------------------------

%after q = 0, which corresponds to jj = 1, we now need to loop from:
jj_max = round(N.Nc/2); 
annuls_smplng = 1;
loop_range = ((1 + annuls_smplng) : annuls_smplng : jj_max);

%-------------------------------------------------------------------------

prtf(2 : (1 + length(loop_range))) = 0;

for jj = loop_range
  
  
  o_diam_circ = (((N.Nc/2+1-xx)./jj).^2 + ((N.Nr/2+1-yy)./jj).^2 < 1);
  i_diam_circ = (((N.Nc/2+1-xx)./(jj-1)).^2 + ((N.Nr/2+1-yy)./(jj-1)).^2 < 1);
  
  temp1 = beam_stop .* (o_diam_circ - i_diam_circ);
  
  prtf(kk) = sum(sum( temp1 .* soln_modulus )) / (eps + sum(sum( temp1 .* meas_modulus )) );

  kk = kk+1;

end

%-------------------------------------------------------------------------

end

