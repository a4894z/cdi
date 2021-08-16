function [ V_bounded ] = phase_limits_project( V, lim )

[ phs_V ] = array_value_limits( angle( V ), lim( 1 ), lim( 2 ) );

V_bounded = abs( V ) .* exp( 1i * phs_V );






