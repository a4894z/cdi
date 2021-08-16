function [ V_bounded_abs ] = scale_matrixcolumns_magnitude( V, absminmax )

abs_V = abs( V );

min_abs = min( abs_V );
abs_V = abs_V - repmat( min_abs, [ size( V, 1 ), 1 ] ); 

max_abs = max( abs_V );
abs_V = abs_V ./ ( 1e-7 + repmat( max_abs, [ size( V, 1 ), 1 ] )); 

abs_Ds_bounded = ( abs_V * ( absminmax( 2 ) - absminmax( 1 )) + absminmax( 1 ));

V_bounded_abs = abs_Ds_bounded .* exp( 1i * angle( V ));
