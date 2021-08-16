function FDxyz = setup_edgedetect_forwarddiff3Dxyz( FDxyz, sz )

%==================================================================================================

% forward difference edge detection convolution kernels

% in the columns direction
FDxyz.Dx = zeros( sz( 1 ), sz( 2 ), sz( 3 ), 'single' ); 
FDxyz.Dx( 1, 1, 1 ) = +1; 
FDxyz.Dx( 1, 2, 1 ) = -1;

% in the rows direction
FDxyz.Dy = zeros( sz( 1 ), sz( 2 ), sz( 3 ), 'single' ); 
FDxyz.Dy( 1, 1, 1 ) = +1; 
FDxyz.Dy( 2, 1, 1 ) = -1;

% in the pages direction
FDxyz.Dz = zeros( sz( 1 ), sz( 2 ), sz( 3 ), 'single' ); 
FDxyz.Dz( 1, 1, 1 ) = +1; 
FDxyz.Dz( 1, 1, 2 ) = -1;

%==================================================================================================

FDxyz.Dx_fft = fft2( FDxyz.Dx );    % the high pass (in Fourier space) filter in the -- (horizontal) x direction
FDxyz.Dy_fft = fft2( FDxyz.Dy );    % the high pass (in Fourier space) filter in the | (vertical) y direction
FDxyz.Dz_fft = fft2( FDxyz.Dz );

%==================================================================================================

% complex conjugates of the hpfs (used in inverse edge detection operator):
FDxyz.Dy_fft_conj = conj( FDxyz.Dy_fft ); 
FDxyz.Dx_fft_conj = conj( FDxyz.Dx_fft ); 
FDxyz.Dz_fft_conj = conj( FDxyz.Dz_fft ); 

%==================================================================================================

% low pass filter to get rid of noise, high spatial frequency goofiness, etc
FDxyz.qLPF = ones( sz( 1 ), sz( 2 ), sz( 3 ), 'single' ); 

% denominator scaling for inverse edge detection operator:
FDxyz.fft_scaling = FDxyz.qLPF ./ ( abs( FDxyz.Dx_fft ).^2 + abs( FDxyz.Dy_fft ).^2 + abs( FDxyz.Dz_fft ).^2 + 1e-7 );

%==================================================================================================

% don't need the edge detection matrices anymore:
FDxyz = rmfield( FDxyz, {'Dx','Dy','Dz'} );

%==================================================================================================
