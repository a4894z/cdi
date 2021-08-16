function [ V ] = lpf_gauss( V, stdev )

%=========================================================================
%------- Low pass filter using Gaussian filter in Fourier space ----------
%=========================================================================

N = size( V );
mu = 0.5 * N + 1;

% stdev = [ stdev( 1 ) * N( 1 ), stdev( 2 ) * N( 2 ) ];

if length( N ) == 2
    
    low_pass = make_2Dgaussian( N, mu, stdev );
    
elseif length( N ) == 3
    
    low_pass = make_3Dgaussian( N, mu, stdev );
    
end

V = ifftn( fftn( V ) .* fftshift( low_pass ));
