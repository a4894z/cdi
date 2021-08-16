function [ V_bounded_phs ] = scale_matrixcolumns_phase( V, phsminmax )

phs_V = angle( V );

min_phs = min( phs_V );
phs_V = phs_V - repmat( min_phs, [ size( V, 1 ), 1 ] ); 

max_phs = max( phs_V );
phs_V = phs_V ./ ( 1e-7 + repmat( max_phs, [ size( V, 1 ), 1 ] ));

phs_V_bounded = phs_V * ( phsminmax( 2 ) - phsminmax( 1 )) + phsminmax( 1 );

V_bounded_phs = abs( V ) .* exp( 1i * phs_V_bounded );