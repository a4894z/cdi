function FDxy = setup_edgedetect_forwarddiff2Dxy( FDxy, sz )

%==================================================================================================

% forward difference edge detection convolution kernels

% in the horizontal direction
FDxy.Dx = zeros( sz, 'single' ); 
% FDxy.Dx(1,1) = +1; 
% FDxy.Dx(1,2) = 0;
% FDxy.Dx(1,3) = -1;
FDxy.Dx(1,1) = +1; 
FDxy.Dx(1,2) = -1;

% in the vertical direction
FDxy.Dy = zeros( sz, 'single' ); 
% FDxy.Dy(1,1) = +1; 
% FDxy.Dy(2,1) = 0; 
% FDxy.Dy(3,1) = -1;
FDxy.Dy(1,1) = +1; 
FDxy.Dy(2,1) = -1;

%==================================================================================================

% low pass filter to get rid of noise, high spatial frequency goofiness, etc
FDxy.qLPF = ones( sz, 'single' ); 

%==================================================================================================

FDxy.Dx_fft = fft2( FDxy.Dx );    % the high pass (in Fourier space) filter in the horizontal direction
FDxy.Dy_fft = fft2( FDxy.Dy );    % the high pass (in Fourier space) filter in the vertical direction

FDxy.sqrt_rc = single( sqrt( prod( sz )));
% FDxy.Dx_fft = fft2( FDxy.Dx ) / FDxy.sqrt_rc;    
% FDxy.Dy_fft = fft2( FDxy.Dy ) / FDxy.sqrt_rc;    

%==================================================================================================

% complex conjugates of the hpfs (used in inverse edge detection operator):
FDxy.Dy_fft_conj = conj( FDxy.Dy_fft ); 
FDxy.Dx_fft_conj = conj( FDxy.Dx_fft ); 

%==================================================================================================

% denominator scaling for inverse edge detection operator:
FDxy.fft_scaling = FDxy.qLPF ./ ( abs( FDxy.Dx_fft ) .^ 2 + abs( FDxy.Dy_fft ) .^ 2 + 1e-7 );

%==================================================================================================

% don't need the edge detection matrices anymore:
FDxy = rmfield( FDxy, { 'Dx','Dy' });

%==================================================================================================
