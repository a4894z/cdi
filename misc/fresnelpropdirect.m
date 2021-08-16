function [ rhof ] = fresnelpropdirect( rhoi, p )

%==================================================================================================

% some defaults, if not defined:

if ~isfield( p, 'extcurv' ), p.extcurv = true; end
if ~isfield( p, 'intcurv' ), p.intcurv = true; end

%==================================================================================================

% 2D coordinate system indices:

sz = size( rhoi );

x2 = ( -0.5 * sz( 2 ) : 0.5 * sz( 2 ) - 1 ) .^ 2; 
y2 = transpose( ( -0.5 * sz( 1 ) : 0.5 * sz( 1 ) - 1 )) .^ 2; 

%==================================================================================================

% quadratic phase curvature at the starting and ending propagation planes

temp1 = +1i * pi / ( p.lambda * p.zif );

if ( logical( p.extcurv ) == true )
    
    A = fftshift( single( exp( y2 * temp1 * p.zf.dr ^ 2 ) * exp( x2 * temp1 * p.zf.dc ^ 2 )));
    
else
    
    A = ones( sz ); 
    
end

if ( logical( p.intcurv ) == true )
    
    B = single( exp( y2 * temp1 * p.zi.dr ^ 2 ) * exp( x2 * temp1 * p.zi.dc ^ 2 ));
    
else
    
    B = ones( sz );
    
end

%==================================================================================================

% propagate from initial plane zi to final plane zf using direct method:

if strcmp( p.dir, 'forward' ) 
    
    rhof = fftshift( A .* fft2( fftshift( B .* rhoi ))) / sqrt( sz( 1 ) * sz( 2 ));
    
elseif strcmp( p.dir, 'backward' ) 
    
	rhof = fftshift( A .* ifft2( fftshift( B .* rhoi ))) * sqrt( sz( 1 ) * sz( 2 )); 
    
end

%==================================================================================================






