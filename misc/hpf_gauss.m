function [ hpf_me, high_pass ] = hpf_gauss( hpf_me, ssigma_x, ssigma_y, xx, yy)

%
% High pass filter using Gaussian filter in Fourier space
%

if ~( exist('xx','var') && exist('yy','var') )

  Sxy = size(hpf_me);
  Ny = Sxy(1); Nx = Sxy(2); 
  [xx,yy] = meshgrid(-Nx/2 : Nx/2-1, -Ny/2 : Ny/2-1);
  xx = single(xx); 
  yy = single(yy); 
  
  clear('Sxy','Nx','Ny')

end

high_pass =  fftshift( exp(-(xx .^ 2)/(2*(eps + ssigma_x)^2)) .* exp(-(yy .^ 2)/(2*(eps + ssigma_y)^2)) );
%high_pass = high_pass / max(max( high_pass ));
high_pass = 1 - high_pass;
high_pass = high_pass / max(max( high_pass ));

%{
high_pass_zero = (high_pass == 0);

high_pass(high_pass_zero) = inf;

min_high_pass = min(high_pass(:));

high_pass(high_pass_zero) = min_high_pass * 0.1;
%}

%high_pass(high_pass == 0) = 1e-3;

%hpf_me = fftshift( ifft2( fft2( fftshift(hpf_me) ) .* high_pass ) );
hpf_me = ifft2( fft2( hpf_me ) .* high_pass );


end

