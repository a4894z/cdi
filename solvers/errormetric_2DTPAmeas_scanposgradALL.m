function [ f, gf ] = errormetric_2DTPAmeas_scanposgrad( P, T, vs, meas, rs, which, shifttype )

f = [];
gf = [];

szT = single( size( T ));
szP  = single( size( P ));
sqrt_rcT = sqrt( szT(1) * szT(2) );
sqrt_rc = sqrt( szP(1) * szP(2) );

% % temp0 = circshift( T, round( rs ));
% temp0 = subpixelshift2D( T, rs );
% phi = temp0( vs.r, vs.c ) .* P;

[ phi, ~ ] = enforce_2DTPAsposview( P, T, vs, rs, shifttype );

%==================================================================================================

Fphi = fft( fftshift( phi, 1 ), [], 1 );
Fphi = fft( fftshift( Fphi, 2 ), [], 2 ) / sqrt_rc;
% Fphi = fft( phi, [], 1 );
% Fphi = fft( Fphi, [], 2 ) / sqrt_rc;
abs_Fphi2 = abs( Fphi ) .^ 2;
sum_abs_Fphi2 = sum( abs_Fphi2, 3 );        % sum over modes

sum_sqrt_abs_Fphi = sqrt( sum_abs_Fphi2 );
% sum_sqrt_abs_Fphi = meas.Dneq0 .* sqrt( sum_abs_Fphi2 );

errmetric_diff = meas.Dneq0 .* ( sum_sqrt_abs_Fphi - meas.D );
% errmetric_diff = meas.Dneq0 .* sum_sqrt_abs_Fphi - meas.D;
% errmetric_diff = ( sum_sqrt_abs_Fphi - meas.D );

% tmp0 = meas.Dneq0 .* ( 1 - meas.D ./ ( 1e-7 + sum_sqrt_abs_Fphi ));

%==================================================================================================

% std measurement metric:

if strcmp( which, 'val_grad' ) || strcmp( which, 'val_gradFD' ) || strcmp( which, 'val' )

    f = sum( abs( errmetric_diff( : )).^2 );
    
end

%==================================================================================================

% grad wrt scan position of std measurement metric:

if strcmp( which, 'val_grad' ) || strcmp( which, 'grad' )
    
    %============================================
    
%     [ qx, qy ] = meshgrid( 1 : szT( 1 ), 1 : 1 : szT( 2 ));
%     qy = fftshift( qy - 0.5 * szT( 1 ) - 1 );
%     qx = fftshift( qx - 0.5 * szT( 2 ) - 1 );

% % tic
% %     [ qx, qy ] = meshgrid( ( -0.5 * szT( 2 ) + 1 ) : ( 0.5 * szT( 2 ) - 0 ), ( -0.5 * szT( 1 ) + 1 ) : ( 0.5 * szT( 1 ) - 0 ));
%     [ qx, qy ] = meshgrid( ( -0.5 * szT( 2 )) : ( 0.5 * szT( 2 ) - 1 ), ( -0.5 * szT( 1 )) : ( 0.5 * szT( 1 ) - 1 ));
%     qx = fftshift( qx );
%     qy = fftshift( qy );
%     
%     
%     
% 
%     phsramp = exp( -2 * pi * 1i * ( qx * rs( 2 ) / szT(2) + qy * rs( 1 ) / szT(1) ));
%     FT = fft2( fftshift( T )) / ( 0 + 1 * sqrt_rcT);
%     Tshft = FT .* phsramp;
% 
% 
%     qcTshft2 = fftshift( ifft2( qx .* Tshft )) * sqrt_rcT;
%     qrTshft2 = fftshift( ifft2( qy .* Tshft )) * sqrt_rcT;
%     
%     
%     qcTshft = qcTshft2( vs.r, vs.c );
%     qrTshft = qrTshft2( vs.r, vs.c );
%     
% % toc
%     
    
    
    

% tic
    X = fft2( T );

    % floors take care of odd-length signals.
    x_shift = exp( -2 * pi * 1i * rs(2) * [ 0 : floor( szT( 2 ) / 2 ) - 1, floor( -szT( 2 ) / 2 ) : -1 ]  / szT( 2 ) );
    y_shift = exp( -2 * pi * 1i * rs(1) * [ 0 : floor( szT( 1 ) / 2 ) - 1, floor( -szT( 1 ) / 2 ) : -1 ]' / szT( 1 ) );

    Tshft = X .* ( y_shift * x_shift );

    qx = repmat( [ 0 : floor( szT( 2 ) / 2 ) - 1, floor( -szT( 2 ) / 2 ) : -1 ], szT( 1 ), 1 );
    qy = repmat( [ 0 : floor( szT( 1 ) / 2 ) - 1, floor( -szT( 1 ) / 2 ) : -1 ].', 1, szT( 2 ));

    qcTshft = ifft2( qx .* Tshft );
    qrTshft = ifft2( qy .* Tshft );

    qcTshft = qcTshft( vs.r, vs.c );
    qrTshft = qrTshft( vs.r, vs.c );
% toc

% qcTshft = qcTshft / max(abs( qcTshft(:)));
% qrTshft = qrTshft / max(abs( qrTshft(:)));

    %============================================


    
    
    
    
%     sum_FPp_iFqxFTFexpxsys = zeros( szP(1), szP(2), 'single');
%     sum_FPp_iFqyFTFexpxsys = zeros( szP(1), szP(2), 'single');
% 
%     for pp = 1 : size( phi, 3 )
% 
%         sum_FPp_iFqxFTFexpxsys = sum_FPp_iFqxFTFexpxsys + ...
%             conj( fft2( fftshift( phi( :, :, pp )))) .* fft2( fftshift( P( :, :, pp ) .* qcTshft )) / sqrt_rcT ^ 2;
%         
%         sum_FPp_iFqyFTFexpxsys = sum_FPp_iFqyFTFexpxsys + ...
%             conj( fft2( fftshift( phi( :, :, pp )))) .* fft2( fftshift( P( :, :, pp ) .* qrTshft )) / sqrt_rcT ^ 2;
% 
%     end

    




%     Fphi = fft( fftshift( phi, 1 ), [], 1 );
%     Fphi = fft( fftshift( Fphi, 2 ), [], 2 );
    
    tmp0x = fft( fftshift( P .* qcTshft, 1 ), [], 1 );
    tmp0x = fft( fftshift( tmp0x, 2 ), [], 2 ) / sqrt(numel(P));
    
    tmp0y = fft( fftshift( P .* qrTshft, 1 ), [], 1 );
    tmp0y = fft( fftshift( tmp0y, 2 ), [], 2 ) / sqrt(numel(P));
    
    sum_FPp_iFqxFTFexpxsys = sum( conj( Fphi ) .* tmp0x, 3 );   
    sum_FPp_iFqyFTFexpxsys = sum( conj( Fphi ) .* tmp0y, 3 );   
   
    
    
    
    
    
 
    tmp0 = errmetric_diff ./ ( 1e-7 + sum_sqrt_abs_Fphi );
%     tmp1 = meas.Dneq0 .* ( 1 - meas.D ./ ( 1e-7 + sum_sqrt_abs_Fphi ));

    tmp0x = tmp0 .* imag( sum_FPp_iFqxFTFexpxsys );
    tmp0y = tmp0 .* imag( sum_FPp_iFqyFTFexpxsys );
    
    gf( 1 ) = ( 4 * pi / szP(1) ) * sum( tmp0y( : ));
    gf( 2 ) = ( 4 * pi / szP(2) ) * sum( tmp0x( : ));
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%     gf = gf / norm( gf );
      
end


%==================================================================================================

% FORWARD FINITE DIFFERENCE GRAD WRT SCAN POS

if strcmp( which, 'val_gradFD' ) || strcmp( which, 'gradFD' )
    
    
    
end

