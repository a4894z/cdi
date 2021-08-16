function [ cf, gf ] = errormetric_2DTPAmeas_scanposgrad_TESTING( P, T, vs, meas, spos, rs0, rsALL, ss, indxsubset, which, shifttype )

cf = [];
gf = [];

% szT = single( size( T ));
% sqrt_rcT = sqrt( szT(1) * szT(2) );

szP  = single( size( P ));
sqrt_rc = sqrt( szP(1) * szP(2) );

% % temp0 = circshift( T, round( spos ));
% temp0 = subpixelshift2D( T, spos );
% phi = temp0( vs.r, vs.c ) .* P;

[ phi, ~ ] = enforce_2DTPAsposview( P, T, vs, spos, shifttype );

%==================================================================================================

% Fphi = fft( fftshift( phi, 1 ), [], 1 ); Fphi = fft( fftshift( Fphi, 2 ), [], 2 ) / sqrt_rc;
Fphi = fft( fft( phi, [], 1 ), [], 2 ) / sqrt_rc;
abs_Fphi2 = abs( Fphi ) .^ 2;
sum_abs_Fphi2 = sum( abs_Fphi2, 3 );        

% sum_sqrt_abs_Fphi = sqrt( sum_abs_Fphi2 );
% errmetric_diff = meas.Dneq0 .* ( sum_sqrt_abs_Fphi - meas.D );
% % errmetric_diff = ( sum_sqrt_abs_Fphi - meas.D );

sum_sqrt_abs_Fphi = sqrt( meas.Dneq0 .* sum_abs_Fphi2 );
errmetric_diff = sum_sqrt_abs_Fphi - meas.D;

% tmp0 = meas.Dneq0 .* ( 1 - meas.D ./ ( 1e-7 + sum_sqrt_abs_Fphi ));

%==================================================================================================

% std measurement metric:

if strcmp( which, 'val_grad' ) || strcmp( which, 'val_gradFD' ) || strcmp( which, 'val' )

    cf.f = sum( abs( errmetric_diff( : )).^2 );
    
%     rsALL = sol.spos.rs;
%     rs0 = expt.spos.rs( ss, : );
%     indxsubset = expt.spos.indxsubset;
    
    % compute 1 / ( rs( aa = ss ) - rs( aa ~= ss )) 
    tmp9 = 0;
    for aa = 1 : length( indxsubset )

        if indxsubset( aa ) ~= indxsubset( ss )

            rsNN = rsALL( aa, : );
            tmp9 = tmp9 + 1 / ( 1e-4 + sum( abs( spos - rsNN ).^1 ));

        end

    end

    cf.fNN = tmp9;

    cf.fMAX = sum( abs( spos - rs0 ).^2 );
    
            
end

%==================================================================================================
