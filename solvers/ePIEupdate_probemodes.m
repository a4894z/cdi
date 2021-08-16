function [ phi ] = ePIEupdate_probemodes( psi, phi, T, vs_r, vs_c, rs, update_order, shifttype )

max_abs_sample_abs2 = max( abs( T( : )) .^ 2 );

for ss = update_order     % order * DOES * matter here !!!

    TFview = getview_2DsampleTF( T, vs_r, vs_c, rs( ss, : ), shifttype );
    
    % ePIE
    phi = phi + conj( TFview ) .* ( psi( :, :, :, ss ) - phi .* TFview ) / max_abs_sample_abs2;
    
    % rPIE
%     aa = 0.25;
%     phi = phi + conj( TFview ) .* ( psi( :, :, :, ss ) - phi .* TFview ) ./ ( aa * max_abs_sample_abs2 + ( 1 - aa ) * abs( TFview ) .^ 2 );
   
end
