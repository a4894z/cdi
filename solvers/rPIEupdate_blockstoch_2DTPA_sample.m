function [ T ] = rPIEupdate_blockstoch_2DTPA_sample( psi,          ...
                                                     phi,          ...
                                                     T,            ...
                                                     vs_r,         ...
                                                     vs_c,         ...
                                                     rs_0,         ...
                                                     update_order, ...
                                                     rPIE_alpha,   ...
                                                     shifttype )

sum_abs_probe_abs2     = sum( abs( phi ) .^ 2, 3 );
max_sum_abs_probe_abs2 = max( sum_abs_probe_abs2( : ));

conj_probe = conj( phi );

for ss = update_order      % order ***DOES*** matter here !!!
  
    rs = +rs_0( ss, : );                 
    
    [ TFview ] = getview_2DsampleTF( T, vs_r, vs_c, rs, shifttype );

    update_term_T = sum( conj_probe .* ( psi( :, :, :, ss ) - phi .* TFview ), 3 );

    update_term_T = update_term_T ./ ( rPIE_alpha * max_sum_abs_probe_abs2 + ( 1 - rPIE_alpha ) * sum_abs_probe_abs2 );
    
    T( round( -1.0 * rs( 1 ) + vs_r ), ...
       round( -1.0 * rs( 2 ) + vs_c ))     = TFview + update_term_T;
      
end

